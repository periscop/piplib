#!/usr/bin/python3
import sys
import os
import subprocess
import string

# Check arguments
if len(sys.argv) < 3:
    print('Usage:', sys.argv[0], 'file.dat file.ll [pip_path]')
    print("  file.dat # input to the parallel integer programming solver")
    print("  file.ll  # expected output from the solver")
    sys.exit(2)

# Get arguments
source_filename = sys.argv[1]
ll_filename = sys.argv[2]
pip = "exemple32" if len(sys.argv) < 4 else sys.argv[3]

# Display arguments
print(sys.argv)

def check_exists(filename):
    if not os.path.exists(filename):
        print(f"Error: file '{filename}' does not exist")
        sys.exit(3)

for f in [source_filename, ll_filename, pip]:
    check_exists(f)

def text_contents(filename):
    with open(filename, "r") as f:
        return f.read()

# Get source
source = text_contents(source_filename)

# Get ll form pip
pip_output = subprocess.Popen([pip, source_filename], shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE)
pip_output_string = pip_output.communicate(input=source.encode())[0].decode()
pip_output.stdin.close()

# Get ll
ll = text_contents(ll_filename)

def whitespace_normalized(s):
    # trim lines and filter out empties
    rv = ""
    for line in s.split('\n'):
        space_normed = ' '.join(line.split())
        if space_normed:
            rv = rv + space_normed
    return rv

# Compare pip_output and ll
s0 = whitespace_normalized(pip_output_string)
s1 = whitespace_normalized(ll)

# Result
if s0 != s1:
    print(f"Result: {pip} < {source_filename} and {ll_filename} are different:")
    print("Expected Output:")
    print(s1)
    print("Actual Output:")
    print(s0)
    sys.exit(1)
else:
    print(f"Result: {pip} < {source_filename} and {ll_filename} are identical")

# End
sys.exit(0)
