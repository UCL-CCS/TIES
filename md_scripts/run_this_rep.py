"""
Run this replica.

For each executed state, update the local file in order to remember the state.
"""

state = open('state').read()

if state == "CREATED":
    with open('state', 'w') as FOUT:
        FOUT.write('EQ_RUNNING')

    # run the NAMD eq
    # how to run NAMD?
    with open('state', 'w') as FOUT:
        FOUT.write('EQ_DONE')

# run the MD simulation