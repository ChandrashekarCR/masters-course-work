import subprocess
import os
import argparse
import ray


def parse_arguments():
    parser = argparse.ArgumentParser(description="Run ADMIXTURE in parallel using Ray.")
    parser.add_argument("--bed_file", required=True, help="Path to the PLINK BED file.")
    parser.add_argument("--log_dir", required=True, help="Directory to store log files.")
    parser.add_argument("--min_k", type=int, required=True, help="Minimum K value.")
    parser.add_argument("--max_k", type=int, required=True, help="Maximum K value.")
    return parser.parse_args()

# Initialize Ray with all available resources
ray.init(num_cpus=64) #num_cpus=os.cpu_count()

@ray.remote
def run_admixture(k, bed_file, log_dir):
    #Runs ADMIXTURE for a given k in parallel.
    log_file = os.path.join(log_dir, f"log{k}.out")
    cmd = f"admixture --cv {bed_file} {k} > {log_file}"
    subprocess.run(cmd, shell=True, check=True)
    return f"Finished ADMIXTURE for K={k}"

if __name__ == "__main__":
    args = parse_arguments()
    
    # Ensure log directory exists
    os.makedirs(args.log_dir, exist_ok=True)
    
    k_values = list(range(args.min_k, args.max_k + 1))
    futures = [run_admixture.remote(k, args.bed_file, args.log_dir) for k in k_values]
    
    # Wait for all tasks to complete
    results = ray.get(futures)
    
    for res in results:
        print(res)
    
    # Shutdown Ray
    ray.shutdown()