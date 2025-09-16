# vdj.py (Modified Version with argparse for Snakemake integration)
import pandas as pd
import numpy as np
from scipy import stats
import warnings
import argparse
import os
warnings.filterwarnings('ignore')

def extract_vdj_usage(file_path, sheets):
    """
    从指定的Excel文件中提取指定sheet的VDJ/VJ组合和频率数据。
    """
    data = []
    for sheet in sheets:
        try:
            df = pd.read_excel(file_path, sheet_name=sheet)
            # 直接使用 VDJ/VJ Combination 和 Frequency 列
            df['Sample'] = sheet
            data.append(df[['VDJ/VJ Combination', 'Frequency', 'Sample']])
        except Exception as e:
            print(f"Warning: Could not read sheet '{sheet}' from {file_path}: {e}")
            continue
    if not data:
        print(f"No data extracted from {file_path}")
        return pd.DataFrame()
    return pd.concat(data, ignore_index=True)

def compute_union_top_vdj_usage(unimm_df, imm_df, top_n=30, min_samples=2, alpha=0.05):
    """
    改进版：计算中位数频率、SD、p-value。
    - min_samples: 只考虑在至少min_samples个样本中出现的VDJ。
    - alpha: 显著性阈值。
    """
    if len(unimm_df) == 0 or len(imm_df) == 0:
        print("Skipping computation due to empty DataFrames")
        return pd.DataFrame()
    
    # 为每个组计算每个样本的VDJ频率（pivot）
    unimm_pivot = unimm_df.pivot_table(index='Sample', columns='VDJ/VJ Combination', values='Frequency', fill_value=0)
    imm_pivot = imm_df.pivot_table(index='Sample', columns='VDJ/VJ Combination', values='Frequency', fill_value=0)
    
    # 获取top_n的并集VDJ
    unimm_top = unimm_pivot.sum().sort_values(ascending=False).head(top_n).index
    imm_top = imm_pivot.sum().sort_values(ascending=False).head(top_n).index
    union_vdjs = set(unimm_top).union(set(imm_top))
    
    data = []
    for vdj in union_vdjs:
        unimm_freqs = unimm_pivot[vdj].values  # 每个样本的频率数组
        imm_freqs = imm_pivot[vdj].values
        
        # 过滤：只如果在min_samples中>0
        if np.sum(unimm_freqs > 0) < min_samples and np.sum(imm_freqs > 0) < min_samples:
            continue
        
        # 中位数（更鲁棒）
        unimm_median = np.median(unimm_freqs)
        imm_median = np.median(imm_freqs)
        
        # 平均和SD（可选，辅助）
        unimm_mean = np.mean(unimm_freqs)
        imm_mean = np.mean(imm_freqs)
        unimm_sd = np.std(unimm_freqs)
        imm_sd = np.std(imm_freqs)
        
        # 差异：用中位数差
        difference = abs(imm_median - unimm_median)
        
        # 统计检验：Mann-Whitney U (非参数，适合小样本频率)
        if len(unimm_freqs) > 1 and len(imm_freqs) > 1 and np.sum(unimm_freqs > 0) > 0 and np.sum(imm_freqs > 0) > 0:
            stat, pvalue = stats.mannwhitneyu(unimm_freqs, imm_freqs, alternative='two-sided')
            significant = pvalue < alpha
        else:
            pvalue = np.nan
            significant = False
        
        data.append({
            'Combination': vdj,
            'Unimm_Median': unimm_median,
            'Unimm_Mean': unimm_mean,
            'Unimm_SD': unimm_sd,
            'Imm_Median': imm_median,
            'Imm_Mean': imm_mean,
            'Imm_SD': imm_sd,
            'Difference_Median': difference,
            'P_value': pvalue,
            'Significant': significant
        })
    
    df = pd.DataFrame(data)
    # 排序：先按中位数差异降序，再按p-value升序
    df = df.sort_values(['Difference_Median', 'P_value'], ascending=[False, True])
    return df

def main():
    parser = argparse.ArgumentParser(description="VDJ Usage Comparison between Unimmunized and Immunized Samples")
    parser.add_argument('--unimm_file', required=True, help="Path to unimmunized top_combinations_summary.xlsx")
    parser.add_argument('--imm_file', required=True, help="Path to immunized top_combinations_summary.xlsx")
    parser.add_argument('--unimm_vh_sheets', nargs='+', required=True, help="List of unimmunized VH sheet names")
    parser.add_argument('--unimm_vk_sheets', nargs='+', required=True, help="List of unimmunized VK sheet names")
    parser.add_argument('--imm_vh_sheets', nargs='+', required=True, help="List of immunized VH sheet names")
    parser.add_argument('--imm_vk_sheets', nargs='+', required=True, help="List of immunized VK sheet names")
    parser.add_argument('--top_n', type=int, default=30, help="Top N combinations for union (default: 30)")
    parser.add_argument('--min_samples', type=int, default=2, help="Minimum samples with >0 frequency (default: 2)")
    parser.add_argument('--alpha', type=float, default=0.05, help="Significance threshold for p-value (default: 0.05)")
    parser.add_argument('--output_dir', default='.', help="Output directory for CSV files (default: current dir)")
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Starting VDJ comparison with config: {vars(args)}")

    # Process VH data
    print("\nProcessing VH data...")
    unimm_vh = extract_vdj_usage(args.unimm_file, args.unimm_vh_sheets)
    imm_vh = extract_vdj_usage(args.imm_file, args.imm_vh_sheets)
    vh_df = compute_union_top_vdj_usage(unimm_vh, imm_vh, args.top_n, args.min_samples, args.alpha)

    # Process VK data
    print("\nProcessing VK data...")
    unimm_vk = extract_vdj_usage(args.unimm_file, args.unimm_vk_sheets)
    imm_vk = extract_vdj_usage(args.imm_file, args.imm_vk_sheets)
    vk_df = compute_union_top_vdj_usage(unimm_vk, imm_vk, args.top_n, args.min_samples, args.alpha)

    if len(vh_df) == 0 and len(vk_df) == 0:
        print("No valid VDJ combinations found for analysis.")
        return

    # Prepare top tables (top 30, but use args.top_n for head)
    if len(vh_df) > 0:
        top_vh = vh_df.head(args.top_n)[['Combination', 'Unimm_Median', 'Imm_Median', 'Difference_Median', 'P_value', 'Significant']]
    else:
        top_vh = pd.DataFrame()

    if len(vk_df) > 0:
        top_vk = vk_df.head(args.top_n)[['Combination', 'Unimm_Median', 'Imm_Median', 'Difference_Median', 'P_value', 'Significant']]
    else:
        top_vk = pd.DataFrame()

    # Save results to CSV files
    vh_top_path = os.path.join(args.output_dir, 'top_vh_differences.csv')
    vk_top_path = os.path.join(args.output_dir, 'top_vk_differences.csv')
    vh_complete_path = os.path.join(args.output_dir, 'complete_vh_analysis.csv')
    vk_complete_path = os.path.join(args.output_dir, 'complete_vk_analysis.csv')

    top_vh.to_csv(vh_top_path, index=False)
    top_vk.to_csv(vk_top_path, index=False)
    vh_df.to_csv(vh_complete_path, index=False)
    vk_df.to_csv(vk_complete_path, index=False)

    # Print results
    print("\nVH: Top VDJ Combinations with Largest Median Differences (Unimmunized -> Immunized)")
    if len(top_vh) > 0:
        print(top_vh.to_string(index=False))
    else:
        print("No VH data available.")

    print("\nVK: Top VDJ Combinations with Largest Median Differences (Unimmunized -> Immunized)")
    if len(top_vk) > 0:
        print(top_vk.to_string(index=False))
    else:
        print("No VK data available.")

    print(f"\nResults saved to:")
    print(f"- {vh_top_path}")
    print(f"- {vk_top_path}")
    print(f"- {vh_complete_path}")
    print(f"- {vk_complete_path}")

if __name__ == "__main__":
    main()
