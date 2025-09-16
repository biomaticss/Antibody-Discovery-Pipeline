# ana_modified.py (Improved with argparse for Snakemake integration)
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from skbio.diversity.alpha import shannon
import gc
import os
import argparse
from matplotlib_venn import venn2
from matplotlib.patches import Patch
import matplotlib as mpl
from collections import defaultdict
from scipy.integrate import trapezoid

# 设置全局字体大小
mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['axes.labelsize'] = 12

def load_mixcr_data(file_path):
    """加载MixCR输出文件并预处理，包含D基因处理"""
    print(f"Loading data from: {file_path}")
    
    usecols = ['cloneId', 'readCount', 'readFraction', 'aaSeqCDR3', 
               'allVHitsWithScore', 'allDHitsWithScore', 'allJHitsWithScore']
    
    dtype_spec = {
        'cloneId': 'str', 
        'readCount': 'int32', 
        'readFraction': 'float32', 
        'aaSeqCDR3': 'str', 
        'allVHitsWithScore': 'str', 
        'allDHitsWithScore': 'str',
        'allJHitsWithScore': 'str'
    }
    
    try:
        chunks = []
        for chunk in pd.read_csv(file_path, sep='\t', usecols=usecols, 
                                dtype=dtype_spec, low_memory=False, chunksize=100000):
            chunks.append(chunk)
            print(f"Loaded chunk with {len(chunk)} rows")
        
        df = pd.concat(chunks, ignore_index=True)
        print(f"Successfully loaded {len(df)} clones from {os.path.basename(file_path)}")
        
        # 提取CDR3长度
        df['CDR3_length'] = df['aaSeqCDR3'].str.len()
        
        # 提取V基因
        df['V_gene'] = df['allVHitsWithScore'].str.split('*').str[0].str.split(',').str[0]
        
        # 提取D基因（如果存在）
        df['D_gene'] = df['allDHitsWithScore'].str.split('*').str[0].str.split(',').str[0]
        df['D_gene'] = df['D_gene'].fillna('')  # 处理缺失值
        
        # 提取J基因
        df['J_gene'] = df['allJHitsWithScore'].str.split('*').str[0].str.split(',').str[0]
        
        # 根据链类型生成组合
        if 'IGH' in file_path:  # 重链使用VDJ组合
            df['VDJ_combination'] = df['V_gene'] + '-' + df['D_gene'] + '-' + df['J_gene']
            df['chain_type'] = 'VH'
        else:  # 轻链使用VJ组合
            df['VDJ_combination'] = df['V_gene'] + '-' + df['J_gene']
            df['chain_type'] = 'VK'
        
        # 计算标准化频率
        total_reads = df['readFraction'].sum()
        if total_reads > 0:
            df['norm_freq'] = df['readFraction'] / total_reads
        else:
            df['norm_freq'] = 0.0
            
        # 删除原始列
        df.drop(columns=['allVHitsWithScore', 'allDHitsWithScore', 'allJHitsWithScore'], inplace=True)
        
        return df
    
    except Exception as e:
        print(f"Error loading file {file_path}: {str(e)}")
        return pd.DataFrame()

def extract_full_length_sequences(df, sample_name, output_dir):
    """提取全长序列并根据readCount复制，保存为纯序列文本文件"""
    if len(df) == 0:
        print(f"Skipping sequence extraction for {sample_name} due to empty DataFrame")
        return
    
    # 选择aaSeqCDR3列并根据readCount复制
    sequences = df.loc[df.index.repeat(df['readCount'])]['aaSeqCDR3']
    
    # 保存到文本文件，不包含列名
    output_file = os.path.join(output_dir, f'full_length_sequences_{sample_name}.txt')
    sequences.to_csv(output_file, index=False, header=False)
    print(f"Saved full-length sequences for {sample_name} to {output_file}")

def get_top_combinations(df, sample_name, n=100, readcount_threshold=3, output_dir='.'):
    """获取并保存样本的前N个VDJ/VJ组合，应用readCount阈值"""
    if len(df) == 0:
        print(f"Skipping top combinations for {sample_name} due to empty DataFrame")
        return None
    
    # 应用readCount阈值
    filtered_df = df[df['readCount'] > readcount_threshold]
    if len(filtered_df) == 0:
        print(f"No clones above readCount threshold for {sample_name}")
        return None
    
    # 按组合分组并计算总频率
    combo_freq = filtered_df.groupby('VDJ_combination')['norm_freq'].sum().reset_index()
    
    # 添加出现次数
    combo_count = filtered_df['VDJ_combination'].value_counts().reset_index()
    combo_count.columns = ['VDJ_combination', 'occurrence_count']
    
    # 合并频率和出现次数
    combo_stats = pd.merge(combo_freq, combo_count, on='VDJ_combination', how='left')
    
    # 添加平均CDR3长度
    cdr3_lengths = filtered_df.groupby('VDJ_combination')['CDR3_length'].mean().reset_index()
    cdr3_lengths.columns = ['VDJ_combination', 'avg_cdr3_length']
    combo_stats = pd.merge(combo_stats, cdr3_lengths, on='VDJ_combination', how='left')
    
    # 排序并获取前N个
    top_combos = combo_stats.sort_values('norm_freq', ascending=False).head(n)
    
    # 添加排名
    top_combos['rank'] = range(1, len(top_combos) + 1)
    
    # 重命名列
    top_combos.columns = ['VDJ/VJ Combination', 'Frequency', 'Occurrence Count', 'Avg CDR3 Length', 'Rank']
    
    # 保存到CSV
    csv_filename = os.path.join(output_dir, f'top_combinations_{sample_name}.csv')
    top_combos.to_csv(csv_filename, index=False)
    print(f"Saved top {n} combinations for {sample_name} to {csv_filename}")
    
    return top_combos

def analyze_replicates(df1, df2, readcount_threshold=3):
    """比较两个样本的重复性（同时计算VDJ组合和CDR3序列的相似性，应用readCount阈值）"""
    if len(df1) == 0 or len(df2) == 0:
        print("Skipping replicate analysis due to empty DataFrame(s)")
        return {}
    
    # 应用readCount阈值
    df1_filtered = df1[df1['readCount'] > readcount_threshold]
    df2_filtered = df2[df2['readCount'] > readcount_threshold]
    
    results = {}
    
    # 1. 基于VDJ组合的相似性分析
    vdj1 = set(df1_filtered['VDJ_combination'])
    vdj2 = set(df2_filtered['VDJ_combination'])
    shared_vdj = vdj1 & vdj2
    
    results['VDJ_Jaccard'] = len(shared_vdj) / len(vdj1 | vdj2) if (vdj1 | vdj2) else 0
    results['VDJ_Shared_ratio1'] = len(shared_vdj) / len(vdj1) if vdj1 else 0
    results['VDJ_Shared_ratio2'] = len(shared_vdj) / len(vdj2) if vdj2 else 0
    
    # 2. 基于CDR3序列的相似性分析
    cdr3_1 = set(df1_filtered['aaSeqCDR3'])
    cdr3_2 = set(df2_filtered['aaSeqCDR3'])
    shared_cdr3 = cdr3_1 & cdr3_2
    
    results['CDR3_Jaccard'] = len(shared_cdr3) / len(cdr3_1 | cdr3_2) if (cdr3_1 | cdr3_2) else 0
    results['CDR3_Shared_ratio1'] = len(shared_cdr3) / len(cdr3_1) if cdr3_1 else 0
    results['CDR3_Shared_ratio2'] = len(shared_cdr3) / len(cdr3_2) if cdr3_2 else 0
    
    # 3. 克隆频率相关性（使用CDR3序列）
    top_df1 = df1_filtered.nlargest(500, 'readCount')
    top_df2 = df2_filtered.nlargest(500, 'readCount')
    
    merged = pd.merge(
        top_df1[['aaSeqCDR3', 'norm_freq']],
        top_df2[['aaSeqCDR3', 'norm_freq']],
        on='aaSeqCDR3',
        how='outer',
        suffixes=('_J1', '_J2')
    ).fillna(0)
    
    if len(merged) > 2:
        spearman_corr = stats.spearmanr(merged['norm_freq_J1'], merged['norm_freq_J2'])
        results['Spearman_corr'] = spearman_corr.correlation
        results['Spearman_pvalue'] = spearman_corr.pvalue
    else:
        results['Spearman_corr'] = np.nan
        results['Spearman_pvalue'] = np.nan
    
    # 可视化
    plt.figure(figsize=(10, 6))
    plt.scatter(
        np.log10(merged['norm_freq_J1'] + 1e-8),
        np.log10(merged['norm_freq_J2'] + 1e-8),
        alpha=0.6
    )
    plt.xlabel('Sample1 Log10(Frequency)')
    plt.ylabel('Sample2 Log10(Frequency)')
    plt.title('Clone Frequency Correlation (Top 500 Clones)')
    plt.savefig(os.path.join(args.output_dir, 'replicate_correlation.png'), dpi=300)  # 使用args.output_dir
    plt.close()
    
    # 释放内存
    del vdj1, vdj2, shared_vdj, cdr3_1, cdr3_2, shared_cdr3, merged
    gc.collect()
    
    return results

def analyze_diversity(df, sample_name, readcount_threshold=3, output_dir='.'):
    """计算抗体多样性指标，应用readCount阈值"""
    if len(df) == 0:
        print(f"Skipping diversity analysis for {sample_name} due to empty DataFrame")
        return {}
    
    # 应用readCount阈值
    filtered_df = df[df['readCount'] > readcount_threshold]
    if len(filtered_df) == 0:
        print(f"No clones above readCount threshold for {sample_name}")
        return {}
    
    results = {}
    
    # 基础统计
    results['Total_clones'] = len(filtered_df)
    results['Unique_CDR3'] = filtered_df['aaSeqCDR3'].nunique()
    results['Unique_VDJ'] = filtered_df['VDJ_combination'].nunique()
    
    # 多样性指数
    if len(filtered_df) > 10000:
        sample_df = filtered_df.sample(10000, random_state=42, weights='readCount')
        results['Shannon_index'] = shannon(sample_df['readCount'])
    else:
        results['Shannon_index'] = shannon(filtered_df['readCount'])
    
    # Gini系数
    if len(filtered_df) > 0:
        sorted_freq = np.sort(filtered_df['norm_freq'])
        cum_wealth = np.cumsum(sorted_freq)
        gini = 1 - 2 * trapezoid(cum_wealth) / (len(cum_wealth) - 1)
        results['Gini_index'] = gini
    else:
        results['Gini_index'] = np.nan
    
    # CDR3长度分布 - 统一横坐标
    plt.figure(figsize=(10, 6))
    sns.histplot(filtered_df['CDR3_length'], bins=30, kde=False)
    plt.title(f'CDR3 Length Distribution - {sample_name}')
    plt.xlim(0, 35)
    plt.savefig(os.path.join(output_dir, f'CDR3_length_{sample_name}.png'), dpi=300)
    plt.close()
    
    # 添加克隆扩增指标
    top10_freq = filtered_df.nlargest(10, 'readCount')['norm_freq'].sum()
    results['Top10_freq'] = top10_freq
    
    return results

def plot_top_genes(df, sample_name, chain_type, readcount_threshold=3, output_dir='.'):
    """绘制每个样本的top10 V基因和J基因使用情况，应用readCount阈值"""
    if len(df) == 0:
        print(f"Skipping top genes plot for {sample_name} due to empty DataFrame")
        return
    
    # 应用readCount阈值
    filtered_df = df[df['readCount'] > readcount_threshold]
    if len(filtered_df) == 0:
        print(f"No clones above readCount threshold for {sample_name}")
        return
    
    # 计算V基因频率
    v_gene_freq = filtered_df.groupby('V_gene')['norm_freq'].sum().reset_index()
    top10_v = v_gene_freq.nlargest(10, 'norm_freq')
    
    # 计算J基因频率
    j_gene_freq = filtered_df.groupby('J_gene')['norm_freq'].sum().reset_index()
    top10_j = j_gene_freq.nlargest(10, 'norm_freq')
    
    # 计算D基因频率（仅对VH链）
    if chain_type == 'VH':
        d_gene_freq = filtered_df.groupby('D_gene')['norm_freq'].sum().reset_index()
        top10_d = d_gene_freq.nlargest(10, 'norm_freq')
    
    # 绘制图形
    plt.figure(figsize=(15, 5) if chain_type == 'VH' else (10, 5))
    
    # 绘制V基因
    ax1 = plt.subplot(1, 3 if chain_type == 'VH' else 2, 1)
    ax1.barh(top10_v['V_gene'], top10_v['norm_freq'], color='skyblue', edgecolor='black')
    ax1.set_title('Top 10 V Genes')
    ax1.set_xlabel('Frequency')
    ax1.set_ylabel('V Gene')
    ax1.invert_yaxis()
    
    # 绘制D基因（仅对VH链）
    if chain_type == 'VH':
        ax2 = plt.subplot(1, 3, 2)
        ax2.barh(top10_d['D_gene'], top10_d['norm_freq'], color='lightcoral', edgecolor='black')
        ax2.set_title('Top 10 D Genes')
        ax2.set_xlabel('Frequency')
        ax2.set_ylabel('D Gene')
        ax2.invert_yaxis()
    
    # 绘制J基因
    ax3 = plt.subplot(1, 3 if chain_type == 'VH' else 2, 2 if chain_type == 'VH' else 2)
    ax3.barh(top10_j['J_gene'], top10_j['norm_freq'], color='lightgreen', edgecolor='black')
    ax3.set_title('Top 10 J Genes')
    ax3.set_xlabel('Frequency')
    ax3.set_ylabel('J Gene')
    ax3.invert_yaxis()
    
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(output_dir, f'top_genes_{sample_name}_{chain_type}.png'), dpi=300)
    plt.close()
    
    print(f"Saved top genes plot for {sample_name}")

def plot_vdj_usage(df, sample_name, chain_type, readcount_threshold=3, output_dir='.'):
    """绘制每个样本的top20 VDJ组合使用情况（水平条形图），应用readCount阈值"""
    if len(df) == 0:
        print(f"Skipping VDJ usage plot for {sample_name} due to empty DataFrame")
        return
    
    # 应用readCount阈值
    filtered_df = df[df['readCount'] > readcount_threshold]
    if len(filtered_df) == 0:
        print(f"No clones above readCount threshold for {sample_name}")
        return
    
    # 计算该样本的VDJ组合频率
    vdj_freq = filtered_df.groupby('VDJ_combination')['norm_freq'].sum().reset_index()
    
    # 取top20组合
    top20 = vdj_freq.nlargest(20, 'norm_freq')
    
    if len(top20) == 0:
        print(f"No VDJ combinations for {sample_name}")
        return
    
    # 按频率排序（从高到低）
    top20 = top20.sort_values('norm_freq', ascending=True)
    
    # 创建水平条形图
    plt.figure(figsize=(12, 10))
    plt.barh(
        top20['VDJ_combination'],
        top20['norm_freq'],
        color='skyblue', 
        edgecolor='black'
    )
    
    plt.title(f'Top 20 VDJ Usage - {sample_name} ({chain_type})')
    plt.xlabel('Frequency')
    plt.ylabel('VDJ Combination')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'VDJ_usage_top20_{sample_name}_{chain_type}.png'), dpi=300)
    plt.close()
    
    # 保存组合信息到CSV
    csv_filename = os.path.join(output_dir, f'vdj_combination_top20_{sample_name}_{chain_type}.csv')
    top20.to_csv(csv_filename, index=False)
    print(f"Saved VDJ combination data for {sample_name} to {csv_filename}")

def plot_group_comparison(df_group1, df_group2, group1_name, group2_name, chain_type, readcount_threshold=3, output_dir='.'):
    """比较两组之间的VDJ组合重叠情况，应用readCount阈值"""
    if len(df_group1) == 0 or len(df_group2) == 0:
        print(f"Skipping group comparison for {group1_name} vs {group2_name} due to empty DataFrame")
        return {}
    
    # 应用readCount阈值
    df_group1_filtered = df_group1[df_group1['readCount'] > readcount_threshold]
    df_group2_filtered = df_group2[df_group2['readCount'] > readcount_threshold]
    
    # 获取唯一的VDJ组合
    vdj_group1 = set(df_group1_filtered['VDJ_combination'])
    vdj_group2 = set(df_group2_filtered['VDJ_combination'])
    
    # 计算重叠统计
    results = {}
    shared = vdj_group1 & vdj_group2
    total_unique = len(vdj_group1 | vdj_group2)
    
    results['Group1_unique'] = len(vdj_group1)
    results['Group2_unique'] = len(vdj_group2)
    results['Shared'] = len(shared)
    results['Jaccard_index'] = len(shared) / total_unique if total_unique > 0 else 0
    
    # 创建韦恩图
    plt.figure(figsize=(8, 6))
    venn2(subsets=(len(vdj_group1 - vdj_group2), 
                   len(vdj_group2 - vdj_group1), 
                   len(shared)),
          set_labels=(group1_name, group2_name))
    plt.title(f'VDJ Combination Overlap - {chain_type}')
    plt.savefig(os.path.join(output_dir, f'venn_{group1_name}_vs_{group2_name}_{chain_type}.png'), dpi=300)
    plt.close()
    
    # 保存统计结果
    stats_df = pd.DataFrame([results])
    stats_df.to_csv(os.path.join(output_dir, f'group_comparison_{group1_name}_vs_{group2_name}_{chain_type}.csv'), index=False)
    
    return results

def main():
    parser = argparse.ArgumentParser(description="Antibody Repertoire Analysis")
    parser.add_argument('--data_dir', required=True, help="MixCR data directory")
    parser.add_argument('--output_dir', required=True, help="Output directory")
    parser.add_argument('--samples', nargs='+', required=True, help="Sample names, e.g., '1175-20-VH 1175-20-VK'")
    parser.add_argument('--readcount_threshold', type=int, default=3, help="Read count threshold")
    parser.add_argument('--top_n', type=int, default=100, help="Top N combinations")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    print("Starting analysis with config:", vars(args))

    # Samples dict (from args)
    samples = {}
    for s in args.samples:
        chain = 'VH' if 'VH' in s else 'VK'
        tsv_file = os.path.join(args.data_dir, f"{s.replace('-VH', '').replace('-VK', '')}_{chain}_assemble_Clones_IG{'H' if chain=='VH' else 'K'}.tsv")
        samples[s] = tsv_file

    # Load data
    data_frames = {}
    for name, path in samples.items():
        print(f"\nLoading {name}...")
        data_frames[name] = load_mixcr_data(path)

    # Extract sequences
    print("\nExtracting full-length sequences...")
    for name, df in data_frames.items():
        extract_full_length_sequences(df, name, args.output_dir)

    # Top combinations
    print("\nExtracting top combinations...")
    top_combinations = {}
    for name, df in data_frames.items():
        top_df = get_top_combinations(df, name, args.top_n, args.readcount_threshold, args.output_dir)
        if top_df is not None:
            top_combinations[name] = top_df

    # Save to Excel (aggregated)
    print("\nSaving all top combinations to Excel...")
    with pd.ExcelWriter(os.path.join(args.output_dir, 'top_combinations_summary.xlsx')) as writer:
        for name, df in top_combinations.items():
            sheet_name = name[:30]  # Excel sheet name limit
            df.to_excel(writer, sheet_name=sheet_name, index=False)
    print(f"Saved to {args.output_dir}/top_combinations_summary.xlsx")

    # Replicate analysis
    print("\nAnalyzing replicates...")
    rep_results = {}
    sample_pairs = [
        ('1175-20-VH', '1175-25-VH'),
        ('1175-20-VK', '1175-25-VK'),
        ('1510-20-VH', '1510-25-VH'),
        ('1510-20-VK', '1510-25-VK')
    ]
    for s1, s2 in sample_pairs:
        if s1 in data_frames and s2 in data_frames:
            print(f"  {s1} vs {s2}...")
            rep_results[f'{s1}_vs_{s2}'] = analyze_replicates(data_frames[s1], data_frames[s2], args.readcount_threshold)
    
    # Save replicate results
    rep_df = pd.DataFrame(rep_results).T
    rep_df.to_csv(os.path.join(args.output_dir, 'replicate_analysis_results.csv'))
    print("\nReplicate Analysis Results:")
    print(rep_df.to_markdown())

    # Diversity analysis
    print("\nAnalyzing diversity...")
    diversity_results = {}
    for name, df in data_frames.items():
        print(f"  Processing {name}...")
        diversity_results[name] = analyze_diversity(df, name, args.readcount_threshold, args.output_dir)
    
    # Save diversity results
    div_df = pd.DataFrame(diversity_results).T
    div_df.to_csv(os.path.join(args.output_dir, 'diversity_metrics.csv'))
    print("\nDiversity Metrics:")
    print(div_df.to_markdown())

    # Plot gene usage and VDJ combinations
    print("\nPlotting gene usage and VDJ combinations...")
    for name, df in data_frames.items():
        print(f"  Plotting {name}...")
        chain_type = 'VH' if 'VH' in name else 'VK'
        plot_top_genes(df, name, chain_type, args.readcount_threshold, args.output_dir)
        plot_vdj_usage(df, name, chain_type, args.readcount_threshold, args.output_dir)

    # Group comparisons
    print("\nAnalyzing group comparisons...")
    comparison_results = []
    group1_vh = pd.concat([data_frames.get('1175-20-VH', pd.DataFrame()), data_frames.get('1175-25-VH', pd.DataFrame())])
    group1_vk = pd.concat([data_frames.get('1175-20-VK', pd.DataFrame()), data_frames.get('1175-25-VK', pd.DataFrame())])
    group2_vh = pd.concat([data_frames.get('1510-20-VH', pd.DataFrame()), data_frames.get('1510-25-VH', pd.DataFrame())])
    group2_vk = pd.concat([data_frames.get('1510-20-VK', pd.DataFrame()), data_frames.get('1510-25-VK', pd.DataFrame())])

    print("  Group1 vs Group2 (VH)...")
    vh_comp = plot_group_comparison(group1_vh, group2_vh, 'Group1', 'Group2', 'VH', args.readcount_threshold, args.output_dir)
    comparison_results.append(vh_comp)

    print("  Group1 vs Group2 (VK)...")
    vk_comp = plot_group_comparison(group1_vk, group2_vk, 'Group1', 'Group2', 'VK', args.readcount_threshold, args.output_dir)
    comparison_results.append(vk_comp)

    print("  Group1 VH vs VK...")
    group1_vh_vk = plot_group_comparison(group1_vh, group1_vk, 'Group1 VH', 'Group1 VK', 'Comparison', args.readcount_threshold, args.output_dir)
    comparison_results.append(group1_vh_vk)

    print("  Group2 VH vs VK...")
    group2_vh_vk = plot_group_comparison(group2_vh, group2_vk, 'Group2 VH', 'Group2 VK', 'Comparison', args.readcount_threshold, args.output_dir)
    comparison_results.append(group2_vh_vk)

    # Save comparison results
    comp_df = pd.DataFrame(comparison_results)
    comp_df.to_csv(os.path.join(args.output_dir, 'group_comparison_results.csv'), index=False)
    print("\nGroup Comparison Results:")
    print(comp_df.to_markdown())

    # Clone expansion visualization
    print("\nPlotting cumulative clone frequency...")
    plt.figure(figsize=(12, 8))
    for name, df in data_frames.items():
        if len(df) > 0:
            top_df = df.nlargest(1000, 'readCount')
            sorted_freq = np.sort(top_df['norm_freq'])[::-1]
            cum_freq = np.cumsum(sorted_freq)
            plt.plot(cum_freq, label=name, linewidth=2)
    
    plt.xlabel('Clone Rank (Top 1000)')
    plt.ylabel('Cumulative Frequency')
    plt.title('Clone Expansion Profile')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, 'cumulative_clone_frequency.png'), dpi=300)
    plt.close()
    
    print("\nAnalysis completed!")

if __name__ == "__main__":
    main()
