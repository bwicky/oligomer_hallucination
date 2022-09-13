# HAL design models

AlphaFold2 models for each sequence that was experimentally tested in the [oligomer hallucination paper](https://www.biorxiv.org/content/10.1101/2022.06.09.493773v1)

For each design, the model with the *highest pTM score* is deposited. All models were generated with 10 recycling steps using `./scoring/AF2.py`.

Each `.pdb` file also contains its AlphaFold2 confidence metrics appended at the end:
- `ptm`: pTM confidence score
- `{plddt|pae|pae_off_diag|pae_contacts_XA}_mean`: mean of the array.
- `{plddt|pae|pae_off_diag|pae_contacts_XA}_median`: median of the array.
- `{plddt|pae|pae_off_diag|pae_contacts_XA}_std`: standard deviation of the array.
- `{plddt|pae|pae_off_diag|pae_contacts_XA}_min`: lowest value of the array.
- `{plddt|pae|pae_off_diag|pae_contacts_XA}_max`: highest value of the array.
- `{plddt|pae|pae_off_diag|pae_contacts_XA}_range`: difference between the max and min of the array.
- `{pae_contacts_XA}_num`: number of contacts.

where `plddt` is the pre-residue pLDDT confidence score returned by AlphaFold2, `pae` is the pAE matrix returned by AlphaFold2, `pae_off_diag` is the subset of *i,j* pairs in the pAE matrix that correspond to inter-chain positions (i.e. ignoring intra-chain positions), and `pae_contacts_XA` is the subset of *i,j* pairs in the pAE matrix that correspond to inter-chain contacts under the angstrom cutoff indicated by `X` (CB-CB distances).
