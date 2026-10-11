"""Run the manuscript verifiers with assertions enabled; requires Python 3.10+ and SymPy."""
import json,subprocess,sys,time
from pathlib import Path
SCRIPTS = ['verify_stratum_symmetry.py', 'verify_global_outcome_classification.py', 'verify_global_linear_reduction.py', 'verify_global_exceptional_reduction.py', 'verify_global_exceptional_complex_term.py', 'verify_global_coarse_reduction.py', 'verify_unified_pair_chart.py', 'verify_coarse_support_bridge.py', 'verify_two_column_kernel.py', 'verify_spectral_moment_rigidity.py', 'verify_triangle_normal_form.py', 'verify_scalar_triangle_exclusion.py', 'verify_commuting_pair_theorem.py', 'verify_fixedpair_structure.py', 'verify_fixedpair_no_extension.py']
def main():
 root=Path(__file__).resolve().parent
 results=[]
 for script in SCRIPTS:
  print(f"Running {script}",flush=True)
  started=time.monotonic()
  try:
   result=subprocess.run([sys.executable,"-u",str(root/script)],cwd=root,capture_output=True,text=True,timeout=180)
   status={"script":script,"exit_code":result.returncode,"seconds":round(time.monotonic()-started,3)}
   if result.returncode:
    print(result.stdout);print(result.stderr,file=sys.stderr)
  except subprocess.TimeoutExpired:
   status={"script":script,"exit_code":None,"timeout_seconds":180}
  results.append(status)
  print(json.dumps(status),flush=True)
  (root/"RUN_RESULTS.json").write_text(json.dumps(results,indent=2)+"\n")
  if status["exit_code"]!=0:return 1
 print("All manuscript verifiers passed.",flush=True)
 return 0
if __name__=="__main__":raise SystemExit(main())
