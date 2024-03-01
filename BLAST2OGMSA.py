# 导入一些python的内置模块和第三方模块
import sys
import argparse
import csv
import statistics

# 定义一个函数，用于读取输入文件的内容，返回一个列表
def read_input_file(input_file):
    data = []
    with open(input_file, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            data.append(row)
    return data

# 定义一个函数，用于处理输入数据，返回一个字典
def process_input_data(data):
    output = {}
    for row in data:
        # 提取各个字段的值
        sample = row[0]
        gene = row[1]
        expression = float(row[2])
        # 根据一些条件进行筛选和计算
        if sample.startswith("Tumor"):
            if gene not in output:
                output[gene] = []
            output[gene].append(expression)
    return output

# 定义一个函数，用于输出结果，打印到标准输出或者指定的输出文件中
def output_result(output, output_file=None):
    # 如果有指定输出文件，就打开它，否则就使用标准输出
    if output_file:
        f = open(output_file, "w")
    else:
        f = sys.stdout
    # 输出表头
    print("Gene\tMean\tMedian\tMin\tMax", file=f)
    # 对输出数据按照基因名字排序
    sorted_output = sorted(output.items(), key=lambda x: x[0])
    # 输出每一行的数据
    for gene, values in sorted_output:
        mean = statistics.mean(values)
        median = statistics.median(values)
        min_value = min(values)
        max_value = max(values)
        print(f"{gene}\t{mean:.2f}\t{median:.2f}\t{min_value:.2f}\t{max_value:.2f}", file=f)
    # 如果有指定输出文件，就关闭它
    if output_file:
        f.close()

# 定义一个函数，用于解析命令行选项，返回一个命名空间对象
def parse_arguments():
    parser = argparse.ArgumentParser(description="ORPA: a tool for analyzing gene expression data")
    parser.add_argument("input_file", help="the path to the input file")
    parser.add_argument("-o", "--output_file", help="the path to the output file (default: stdout)")
    parser.add_argument("-v", "--version", action="version", version="ORPA 1.0")
    args = parser.parse_args()
    return args

# 主程序
def main():
    # 解析命令行选项
    args = parse_arguments()
    # 读取输入文件
    data = read_input_file(args.input_file)
    # 处理输入数据
    output = process_input_data(data)
    # 输出结果
    output_result(output, args.output_file)

# 如果是直接运行这个脚本，就调用主程序
if __name__ == "__main__":
    main()
