import json, os, sys
sys.path.insert(0, "src")
from recurscan.data.dataset import assemble
from recurscan.benchmark import run_one_seed
ds = assemble()
for model in ("coxgcn", "coxcnn"):
    path = f"results/{model}_seeds.jsonl"
    done = set()
    if os.path.exists(path):
        done = {json.loads(l)["seed"] for l in open(path)}
    for s in range(5):
        if s in done:
            continue
        m = run_one_seed(ds, model, s)
        with open(path, "a") as fh:
            fh.write(json.dumps({"seed": s, **m}) + "\n")
        print(model, "seed", s, "C", round(m["c_index"], 4), flush=True)
print("ALL DONE", flush=True)
