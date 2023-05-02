"""
Class that does the wrapper for DSSR execution
"""
import argparse
import os
import subprocess
from typing import Dict, List, Optional, Union
import pandas as pd

from src.utils import save_json


class DSSRWrapper:
    def __init__(
        self,
        dssr_analyse_bin_path: str = os.path.join("dssr", "bin", "analyze"),
        *args,
        **kwargs,
    ):
        """
        :param dssr_analyse_bin_path: path to the dssr analyse binary file
        """
        self.dssr_analyse_bin_path = dssr_analyse_bin_path
        self.setup_env_variables()

    def setup_env_variables(self):
        """
        Setup the X3DNA variable to execute the DSSR code
        """
        base_name = os.path.dirname(os.path.dirname(self.dssr_analyse_bin_path))
        os.environ['X3DNA'] = os.path.join(os.getcwd(), base_name)
        os.environ['PATH'] = f"{os.environ['X3DNA']}/bin:{os.environ['PATH']}"

    def _convert_path_to_list_pdb_files(self, input_path: str) -> List:
        """
        Return a list of all the .pdb files in the folder
        :param input_path: path to a .pdb file or a folder
        """
        if input_path.endswith(".pdb"):
            return [input_path]
        else:
            return [
                os.path.join(input_path, file_)
                for file_ in os.listdir(input_path)
                if file_.endswith(".pdb")
            ]

    def run(
        self, input_path: str, output_path: Optional[str] = None, to_csv: bool = False
    ) -> Dict:
        """
        Get all the angles for a .pdb file (torsion and pseudo-torsion)
        :param input_path: path to a .pdb file or a folder
        :param output_path: a path where to save the files
        :param to_csv: if True, save the output in a .csv file. Otherwise, save it in a .json file
        """
        if input_path is None:
            raise ValueError("input_path must be specified")
        files = self._convert_path_to_list_pdb_files(input_path)
        output = self.get_angles_from_dssr(files, to_csv)
        if output_path is not None:
            if to_csv:
                output.to_csv(output_path, index=False)
            else:
                save_json(output, output_path)
        return output

    def get_angles_from_dssr(self, files: List[str], to_csv: bool) -> Dict:
        """
        Get all the angles for a .pdb file (torsion and pseudo-torsion)
        :param files: list of .pdb files
        :param to_csv: if True, return a dataframe. Otherwise, return a dictionary.
        """
        output = pd.DataFrame({}) if to_csv else {}
        for file_ in files:
            name = os.path.basename(file_).replace(".pdb", "")
            c_output = self.get_all_angles_from_dssr_one_file(file_, to_csv)
            if to_csv:
                output = self._concatenate_dfs([output, c_output], name)
            else:
                output[name] = c_output
        return output

    def _concatenate_dfs(self, dfs: List[pd.DataFrame], name: str) -> pd.DataFrame:
        """
        Concatenate a list of dataframes
        :param dfs: list of dataframes
        :param name: name of the last .pdb file
        :return: a dataframe
        """
        dfs[-1]["name"] = name
        df = pd.concat(dfs)
        df.reset_index(drop=True, inplace=True)
        return df

    def _convert_pd_to_dict(self, df: pd.DataFrame) -> Dict:
        """
        Convert a dataframe to a dictionary
        :param df: a dataframe, output of the C DSSR code
        :return: dictionary of the form
                { 'sequence' : "AC...",
                   'angles': {
                      'eta': [163, 58, ...],
                      'epsilon' : [177, 23, ...],
                      ...
                    }
                }
        """
        columns = [
            "alpha",
            "beta",
            "beta",
            "gamma",
            "delta",
            "epsilon",
            "zeta",
            "eta",
            "theta",
        ]
        angles = {}
        for column in columns:
            angles[column] = df[column].tolist()
        df["sequence"] = df["base"].apply(lambda x: x[-1])
        sequence = "".join(df["sequence"].tolist())
        output = {"sequence": sequence, "angles": angles}
        return output

    def get_all_angles_from_dssr_one_file(
        self, file_: str, to_csv: bool
    ) -> Union[Dict, pd.DataFrame]:
        """
        Get all the torsion and pseudo torsion angles for a .pdb file
        :param file_: path to a .pdb file
        :param to_csv: if True, return a dataframe. Otherwise, return a dictionary.
        :return: a dictionary with the sequence and angles.
        """
        df_torsion = self.get_angles_from_dssr_one_file(file_, is_pseudo=False)
        df_pseudo_torsion = self.get_angles_from_dssr_one_file(file_, is_pseudo=True)
        df_pseudo_torsion.drop(columns=["rank", "base"], inplace=True)
        df = pd.concat((df_torsion, df_pseudo_torsion), axis=1)
        # Add the sequence
        df["sequence"] = df["base"].apply(lambda x: x[-1])
        if to_csv:
            return df
        return self._convert_pd_to_dict(df)

    def get_angles_from_dssr_one_file(
        self, file_: str, is_pseudo: bool = False
    ) -> pd.DataFrame:
        """
        Get all the angles for a .pdb file (torsion and pseudo-torsion)
        :param file_: path to a .pdb file
        :param is_pseudo: whether to output pseudo-torsion angles or not
        """
        c_arg = "-p" if is_pseudo else "-t"
        command = f"{self.dssr_analyse_bin_path} {c_arg} {file_}"
        output = subprocess.check_output(command, shell=True)
        input_string = [line.split(",") for line in output.decode().split("\n")][:-1]
        df = pd.DataFrame(input_string[1:], columns=input_string[0])
        return df

    @staticmethod
    def get_arguments():
        parser = argparse.ArgumentParser(description="DSSR wrapper")
        parser.add_argument(
            "--input_path", type=str, help="path to a .pdb file or a folder"
        )
        parser.add_argument(
            "--output_path", type=str, help="path to save the output in .json or .csv"
        )
        parser.add_argument(
            "--dssr_analyse_bin_path",
            type=str,
            default=os.path.join("dssr", "bin", "analyze"),
            help="path to the dssr analyse binary file",
        )
        parser.add_argument(
            "--to_csv",
            type=bool,
            action=argparse.BooleanOptionalAction,
            help="whether to return the output in a .csv file or not",
        )
        return parser.parse_args()


if __name__ == "__main__":
    args = DSSRWrapper.get_arguments()
    dssr_wrapper = DSSRWrapper(**vars(args))
    dssr_wrapper.run(args.input_path, args.output_path, args.to_csv)
