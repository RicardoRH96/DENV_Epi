import itertools
import subprocess
import sys
import os

def run_collect_times(burn_in, from_state, to_state, input_file):
    """Function to run the collect_times.py script for a given pair of states."""
    # Build the command for each pair of states
    cmd = [
        'python3', 'collect_times.py', str(burn_in), from_state, to_state
    ]
    
    # Open the input file and pass it through stdin to the script
    with open(input_file, 'r') as infile:
        result = subprocess.run(cmd, stdin=infile, capture_output=True, text=True)
        
    # Return the captured stdout output (i.e., the results of the script)
    if result.returncode == 0:
        return result.stdout
    else:
        print(f"Error running collect_times.py for {from_state} -> {to_state}", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        return None

def main(burn_in, states_list, input_file, output_file):
    """Main function to run pairwise transitions for all states."""
    
    # Open the output file to write all transitions
    with open(output_file, 'w') as outfile:
        outfile.write("state\tfrom\tto\ttime\n")  # Write header to the file
        
        # Generate all pairwise combinations of states
        for from_state, to_state in itertools.permutations(states_list, 2):
            print(f"Processing transition: {from_state} -> {to_state}")
            
            # Run the collect_times script for the current pair
            output = run_collect_times(burn_in, from_state, to_state, input_file)
            
            if output:
                # Append the output to the file
                outfile.write(output)

if __name__ == "__main__":
    # Example usage:
    # python3 pairwise_wrapper.py <burn_in> <input_file> <output_file> state1 state2 state3 ...
    
    if len(sys.argv) < 5:
        print("Usage: pairwise_wrapper.py <burn_in> <input_file> <output_file> <state1> <state2> ...")
        sys.exit(-1)
    
    # Parse command-line arguments
    burn_in = int(sys.argv[1])
    input_file = sys.argv[2]
    output_file = sys.argv[3]
    states_list = sys.argv[4:]  # The rest are the states

    # Ensure the collect_times script exists in the same directory
    if not os.path.exists("collect_times.py"):
        print("Error: collect_times.py script not found in the current directory.")
        sys.exit(-1)

    # Call the main function
    main(burn_in, states_list, input_file, output_file)