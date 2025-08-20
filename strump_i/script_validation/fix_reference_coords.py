#!env /usr/bin/python3

import os
import pandas as pd

os.getcwd()

# test parsing.
def parse_line(line):
    linestart = line[:30]
    coord_x = float(line[30:38].strip())
    coord_y = float(line[38:46].strip())
    coord_z = float(line[46:54].strip())
    lineend = line[54:]
    returndict = {"linestart": linestart, "coord_x": coord_x, "coord_y": coord_y, "coord_z": coord_z,
                  "lineend": lineend}
    return returndict

def change_coord(df):
    if df["coord_x"].min() <= -100:
        df["coord_x"] += 100
    if df["coord_y"].min() <= -100:
        df["coord_y"] += 100
    if df["coord_z"].min() <= -100:
        df["coord_z"] += 100
    # float decimal matching.
    df["coord_x"] = [f"{x:.3f}" for x in df["coord_x"]]
    df["coord_y"] = [f"{x:.3f}" for x in df["coord_y"]]
    df["coord_z"] = [f"{x:.3f}" for x in df["coord_z"]]
    return df

def write_line(rows, f):
    writeline = f"{rows[1]}{rows[2]:>8}{rows[3]:>8}{rows[4]:>8}{rows[5]}"
    f.write(writeline)


def parse_pdb(input_path, output_path):
    with open(input_path) as f:
        pdb_df = [parse_line(line) for line in f.readlines()]
    pdb_df = pd.DataFrame(pdb_df)
    pdb_df = change_coord(pdb_df)
    with open(output_path, 'w') as f:
        [write_line(rows, f) for rows in pdb_df.itertuples()]
    pass




def check_pdb(input_path):
    with open(input_path) as f:
        pdb_df = [parse_line(line) for line in f.readlines()]
    pdb_df = pd.DataFrame(pdb_df)
    if pdb_df["coord_x"].min() <= -100 or pdb_df["coord_y"].min() <= -100 or pdb_df["coord_z"].min() <= -100:
        print(os.path.basename(input_path))




pdb_ref_dir = "./reference/renumbered_pdb"
pdb_out_dir = "./reference/pdb"
tgt_list = os.listdir(pdb_ref_dir)

for tgt in tgt_list:
    parse_pdb(os.path.join(pdb_ref_dir, tgt), os.path.join(pdb_out_dir, tgt))



for tgt in tgt_list:
    check_pdb(os.path.join(pdb_ref_dir, tgt))

for tgt in tgt_list:
    check_pdb(os.path.join(pdb_out_dir, tgt))


test = 9.24
test2 = f"{test:.3f}"
test2
f"{test2:>8}"



input_path = os.path.join(pdb_out_dir, tgt)



