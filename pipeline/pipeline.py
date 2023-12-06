# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Functions for building the input features for the AlphaFold model."""

import os
from typing import Any, Mapping, Optional
from absl import logging
from . import parsers
from .jackhmmer import Jackhmmer


def run_msa_tool(
    msa_runner,
    input_fasta_path: str,
    msa_out_path: str,
    msa_format: str,
    use_precomputed_msas: bool,
    max_sto_sequences: Optional[int] = None,
) -> Mapping[str, Any]:
    """Runs an MSA tool, checking if output already exists first."""
    if not use_precomputed_msas or not os.path.exists(msa_out_path):
        if msa_format == "sto" and max_sto_sequences is not None:
            result = msa_runner.query(input_fasta_path, max_sto_sequences)[
                0
            ]  # pytype: disable=wrong-arg-count
        else:
            result = msa_runner.query(input_fasta_path)[0]
        with open(msa_out_path, "w") as f:
            f.write(result[msa_format])
    else:
        logging.warning("Reading MSA from file %s", msa_out_path)
        if msa_format == "sto" and max_sto_sequences is not None:
            precomputed_msa = parsers.truncate_stockholm_msa(
                msa_out_path, max_sto_sequences
            )
            result = {"sto": precomputed_msa}
        else:
            with open(msa_out_path, "r") as f:
                result = {msa_format: f.read()}
    return result


class DataPipeline:
    """Runs the alignment tools and assembles the input features."""

    def __init__(
        self,
        jackhmmer_binary_path: str,
        small_bfd_database_path: Optional[str],
        use_small_bfd: bool,
        use_precomputed_msas: bool = False,
    ):
        """Initializes the data pipeline."""
        self._use_small_bfd = use_small_bfd
        self.jackhmmer_small_bfd_runner = Jackhmmer(
            binary_path=jackhmmer_binary_path,
            database_path=small_bfd_database_path
        )

        self.use_precomputed_msas = use_precomputed_msas

    def process(self, input_fasta_path: str, msa_output_dir: str):
        """Runs alignment tools on the input sequence and creates features."""

        bfd_out_path = os.path.join(msa_output_dir, "small_bfd_hits.sto")
        jackhmmer_small_bfd_result = run_msa_tool(
            msa_runner=self.jackhmmer_small_bfd_runner,
            input_fasta_path=input_fasta_path,
            msa_out_path=bfd_out_path,
            msa_format="sto",
            use_precomputed_msas=self.use_precomputed_msas,
        )
        bfd_msa = parsers.parse_stockholm(jackhmmer_small_bfd_result["sto"])

        return bfd_msa
