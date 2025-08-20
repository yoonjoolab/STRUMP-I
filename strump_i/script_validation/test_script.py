
import argparse
from pathlib import Path
Path(__file__).absolute()

Path(__file__).joinpath("test", "test2")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog='STRUMP-I')
    parser.add_argument("--data_dir", type=str)
    args = parser.parse_args()
    print(args.data_dir)

