import sys

# Get the number of arguments passed
num_args = len(sys.argv)
if num_args not in [2, 4]:
    print("Usage: collect_times <burning #> [from_state] [to_state]")
    sys.exit(-1)

# Burn-in period
burn_in = int(sys.argv[1])

# If three arguments are provided
if num_args == 4:
    print_all = False
    from_state = sys.argv[2]
    to_state = sys.argv[3]
else:
    print_all = True
    from_state = "ALL"
    to_state = "ALL"

# Print diagnostic information to standard error
print(f"Removing states < {burn_in} as burn in.", file=sys.stderr)
print(f"Recording jumps from {from_state} to {to_state}.", file=sys.stderr)

# Skip first 3 lines of input
for _ in range(3):
    input()

state_count = 0
jump_count = 0

# Print header
print("state\tfrom\tto\ttime")

# Process the input line by line
for line in sys.stdin:
    line = line.strip()
    parts = line.split('\t', 2)  # Split the line into three parts
    state = int(parts[0])
    if state >= burn_in:
        state_count += 1
        jumpf = parts[2]
        
        # Remove redundant total count at the end of the jump string
        jumpf = jumpf.rsplit(' ', 1)[0]
        
        # Split jumps and process each jump
        jumps = jumpf.split('},{')
        for jump in jumps:
            jump = jump.replace('{{', '').replace('}}', '')
            site, time, ori, dest = jump.split(',')

            if print_all or (ori == from_state and dest == to_state):
                jump_count += 1
                print(f"{state}\t{ori}\t{dest}\t{time}")

# Print diagnostic information to standard error
print(f"total state counts = {state_count}, total jumps = {jump_count}", file=sys.stderr)