"""查看原始研究资料摘要，或执行条件最大扫描并落盘结果。"""
import argparse
import json
from .io import load
from .domain import summarize
from .scan import analyze, save_result, ExactBudgetExceeded


def main():
    parser = argparse.ArgumentParser(description='合成研究资料工具')
    sub = parser.add_subparsers(dest='command', required=True)

    inspect_parser = sub.add_parser('inspect', help='只显示已有摘要，不修改来源')
    inspect_parser.add_argument('file')

    scan_parser = sub.add_parser('scan', help='条件最大扫描分析（候选族固定在先）')
    scan_parser.add_argument('file')
    scan_parser.add_argument('--outdir', required=True,
                             help='保存候选族、配置摘要、null 最大值与观察排名的目录')
    scan_parser.add_argument('--mode', choices=['auto', 'monte-carlo', 'exact'],
                             default='auto')
    scan_parser.add_argument('--engine', choices=['python', 'numpy'],
                             default='python',
                             help='蒙特卡洛引擎；numpy 为可选加速，候选族与分布与基准一致')

    args = parser.parse_args()
    value, digest = load(args.file)
    if args.command == 'inspect':
        print(json.dumps({'source_sha256': digest, 'summary': summarize(value)},
                         ensure_ascii=False, indent=2, allow_nan=False))
        return

    result = analyze(value, mode=args.mode, engine=args.engine)
    result['source'] = {'file': args.file, 'sha256': digest}
    written = save_result(result, args.outdir)
    print(json.dumps({
        'status': result['status'],
        'mode': result['mode'],
        'family_size': result['family_size'],
        'written': {key: str(path) for key, path in written.items()},
        'notes': result['notes'],
    }, ensure_ascii=False, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
