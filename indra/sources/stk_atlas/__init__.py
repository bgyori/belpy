"""This module implements and API and processor for the atlas of
substrate specificities for the human serine/threonine kinome:

Johnson, J.L., Yaron, T.M., Huntsman, E.M. et al.
An atlas of substrate specificities for the human serine/threonine kinome.
Nature (2023). https://doi.org/10.1038/s41586-022-05575-3

Two files are needed to process the result into INDRA Statements.
- Supplementary Table 3 contains annotations for the Ser/Thr phosphoproteome.
This is a 370 MB xlsx file with rows corresponding to phosphosites and
columns to kinases. Due to issues with "Strict OOXML" processing using
openpyxl/pandas, this file cannot directly be read in Python. To prepare this
file for processing, manually open it with Microsoft Excel, then save the
worksheet "Supplementary Table 3" as a csv file.
- Supplementary Table 1 contains information on each kinase, and due to
non-standard naming in the annotations file, this is needed to process
kinases correctly. This is a small ~43KB xlsx file that can be read directly
by openpyxl/pandas from the web or locally.

For each phosphosite, kinases are ranked 1-303 based on how specific their
substrate motif is for the peptide corresponding to the phosphosite. This
means that one has to choose a rank threshold below which a kinase is
considered to be targeting a given phosphosite. In Johnson et al., for
enrichment analysis the top 15 kinases are taken for each substrate.

"""

from .api import *
