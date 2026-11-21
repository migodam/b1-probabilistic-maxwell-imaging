"""Look up an existing study and its saved input. This command never generates data."""
import argparse
import csv
import os
from pathlib import Path


HERE = Path(__file__).resolve().parent


def frozen_study(root, round_number):
    """Find an archived MAX-STRUCT round in a full B1 checkout or Inv_SLAM."""
    study = Path("research") / f"research_max_struct_r{round_number}"
    public = root / study
    if public.is_dir():
        return public
    return root / "research/delegated/max_struct_r1" / f"r{round_number}_worktree" / study


def read_rows(file):
    with (HERE / file).open(encoding="utf8", newline="") as stream:
        return list(csv.DictReader(stream))


def find_input(row, relpath, workspace, repo):
    bases = []
    if workspace:
        bases.append(workspace / row["workspace_dir"])
    if repo and row["repo_dir"]:
        bases.append(repo / row["repo_dir"])
    return next((base / relpath for base in bases if (base / relpath).is_file()), None)


def main():
    parser = argparse.ArgumentParser(description="查冻结数据，默认列出全部数据族。")
    parser.add_argument("--operator", metavar="ID", help="定位已有算子入口；以 list 列出全部算子。")
    parser.add_argument("dataset", nargs="?")
    parser.add_argument("observation", nargs="?", help="完整观测 ID（包括视图）。")
    parser.add_argument("--list-views", action="store_true", help="列观测 ID、材料聚类键及版本。")
    parser.add_argument("--workspace", type=Path, default=os.environ.get("INV_SLAM_ROOT"))
    parser.add_argument("--repo", type=Path, default=os.environ.get("B1_ROOT"))
    args = parser.parse_args()
    workspace = args.workspace or next((p for p in HERE.parents if (p / "Gaussian").is_dir()), None)
    repo = args.repo or HERE.parents[1]
    if args.operator:
        if args.dataset or args.observation or args.list_views:
            parser.error("算子查询与数据查询分开使用")
        operators = {r["operator"]: r for r in read_rows("OPERATORS.csv")}
        if args.operator == "list":
            for key, op in operators.items():
                print(f'{key:22} {op["device"]} · {op["entry"]}')
            return 0
        if args.operator not in operators:
            parser.error("未知算子；用 --operator list 查看")
        op = operators[args.operator]
        print(f'{op["entry"]}   {op["device"]}\n适用　{op["workload"]}\n费用　{op["accounting"]}')
        candidates = [root / op[field] for root, field in
                      [(workspace, "workspace_file"), (repo, "repo_file")] if root and op[field]]
        available = next((p for p in candidates if p.is_file()), None)
        if available:
            print(available)
            return 0
        print(f'本机缺少源码；B1 仓库位置　{op["repo_file"]}')
        return 1
    datasets = {r["dataset"]: r for r in read_rows("COHORTS.csv")}
    if not args.dataset:
        if args.list_views:
            parser.error("--list-views 需要数据族")
        for name, row in datasets.items():
            print(f'{name:26} {row["origin"]} · {row["population"]}')
        return 0
    if args.dataset not in datasets:
        parser.error("未登记的数据族。无参数运行可列出全部 ID。")
    row = datasets[args.dataset]
    views = [r for r in read_rows("OBSERVATIONS.csv") if r["dataset"] == args.dataset]
    if args.observation:
        view = next((r for r in views if r["observation_id"] == args.observation), None)
        if view is None:
            parser.error("此观测未登记在指定数据族。使用 --list-views 查看版本。")
        path = find_input(row, view["path"], workspace, repo)
        print(f'{view["dataset"]} / {view["observation_id"]}   {view["version"]}')
        print(f'用途　{row["role"]}\n输入　{row["online_fields"]}')
        if path:
            print(path)
            return 0
        print(f'本机缺少 {view["path"]}\n获取　{row["access"]}\n来源　{row["source"]}')
        return 1
    if args.list_views:
        for v in views:
            print(f'{v["observation_id"]}\t{v["object_id"] or "聚类待核"}\t{v["version"]}')
        return 0
    for label, key in [("来源","origin"),("对象","population"),("物理","contract"),
                       ("输入","online_fields"),("用途","role"),("获取","access"),("依据","source")]:
        print(f'{label}　{row[key]}')
    reference = find_input(row, row["reference"], workspace, repo)
    if reference:
        print(f'本地合同　{reference}')
    if row["offline_scoring"]:
        print(f'离线评分代码或资料（相对研究目录）　{row["offline_scoring"]}')
    print(f'登记视图　{len(views)}')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
