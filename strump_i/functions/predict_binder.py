
'''
    STRUMP-I prediction scripts.


'''

# import library
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path

# for make model.
import lightgbm as lgbm

# for metric
from sklearn.metrics import accuracy_score, precision_recall_curve, roc_auc_score, roc_curve, auc

# import manual scaler for data scaling.
from strump_i.functions.strumpi_scaler import StrumpiScaler


# define class
class StrumpiBindingClassifier:
    def __init__(self, model_dir: str):
        self.info_cols = ['allele', 'peptide']

        self.model_path = Path(model_dir, 'lgbm_model_all.pkl')
        with open(self.model_path, 'rb') as f:
            self.model = pickle.load(f)
        self.scaler_path = Path(model_dir, 'lgbm_scaler_all.pkl')
        with open(self.scaler_path, 'rb') as f:
            self.scaler = pickle.load(f)
        self.scaler.startcol=len(self.info_cols)
        self.scaling_features = list(self.scaler.scale_dict.keys())
        self.model_features = self.model.feature_name_
        pass

    def load_data(self, input_dir: str):
        input_path = Path(input_dir)
        quality_path = input_path.joinpath('quality_data.csv')
        energy_path = input_path.joinpath('peptide_energy_data.csv')
        quality_df = pd.read_csv(quality_path)
        energy_df = pd.read_csv(energy_path)
        return_df = pd.merge(energy_df, quality_df, on=['Pdb', 'allele', 'peptide'])
        return_df =  return_df[return_df['final_e'] <= 0]
        return return_df

    def make_pred_df_classifier(self, input_df: pd.DataFrame) -> pd.DataFrame:
        # for classification model.
        pred = self.model.predict_proba(input_df[self.model_features])
        pred = pd.DataFrame(pred)
        pred_df = pd.DataFrame(
            data={
                "allele": input_df['allele'].tolist(),
                'peptide': input_df['peptide'].tolist(),
                'pred': pred[1].tolist()
            }
        )
        return pred_df

    def __call__(self, input_dir, output_dir):
        data = self.load_data(input_dir)
        data = data[self.info_cols + self.scaling_features]
        scaled_data = self.scaler.scaling(data)
        scaled_data.columns = [x.replace(" ", "_") for x in scaled_data.columns]
        data_pred = self.make_pred_df_classifier(scaled_data)
        data_pred.to_csv(Path(output_dir, "pred_result.csv"), header=True, index=False)
        pass


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_dir")
    parser.add_argument("--output_dir")
    parser.add_argument("--model_dir", default=Path(__file__).parent.parent.joinpath("reference", "model", "binding"))
    args = parser.parse_args()
    print(args.model_dir)
    runner = StrumpiBindingClassifier(args.model_dir)
    runner(args.input_dir, args.output_dir)


