conda env config vars set PYTHONNOUSERSITE=1 -n 094   # 094 ignores ~/.local from now on
conda activate 094
aws --version                                          # should say aws-cli/2.38.3
aws configure --profile seqdata                        # paste the two keys, region eu-central-1, output json