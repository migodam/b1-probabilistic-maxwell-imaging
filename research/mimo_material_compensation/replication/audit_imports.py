"""Inventory direct imports and file-read expressions for each study script."""
from pathlib import Path
import ast
import json

HERE=Path(__file__).resolve().parent
STUDY=HERE.parent
GROUPS=("controls","paths","budgets","environment","topology","inputs","replication")
HISTORICAL={"model","shared","verified_wave","native_replay","aperture_control","aperture_joint"}
FILES=sorted(p for group in GROUPS for p in (STUDY/group).rglob("*.py") if "__pycache__" not in p.parts)
LOCAL={p.stem for p in FILES}


def inventory(path):
    tree=ast.parse(path.read_text(),filename=str(path))
    imports=set()
    reads=[]
    for node in ast.walk(tree):
        if isinstance(node,ast.Import):
            imports.update(q.name for q in node.names)
        elif isinstance(node,ast.ImportFrom) and node.module:
            imports.add(node.module)
        elif isinstance(node,ast.Call):
            op=ast.unparse(node.func)
            tail=op.rsplit(".",1)[-1]
            if tail in ("load","read_text","read_bytes"):
                reads.append(dict(line=node.lineno,function=op,
                    argument=ast.unparse(node.args[0]) if node.args else ast.unparse(node.func.value) if isinstance(node.func,ast.Attribute) else ""))
    return dict(path=str(path.relative_to(STUDY)),imports=sorted(imports),
                prior_imports=sorted(m for m in imports if m.split(".")[0] in HISTORICAL),
                project_imports=sorted(m for m in imports if m.split(".")[0] in LOCAL),
                input_expressions=sorted(reads,key=lambda r:r["line"]))


if __name__=="__main__":
    result=[inventory(p) for p in FILES]
    (HERE/"script_dependencies.json").write_text(json.dumps(dict(files=len(result),scripts=result),indent=2,ensure_ascii=False)+"\n")
    print("Inventoried",len(result),"Python scripts")
