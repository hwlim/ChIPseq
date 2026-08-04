#!/usr/bin/env python3
"""Write findDiffPeak.DESeq2.yml for group-wise peak calling.

Given a Group name and sample.tsv, assemble the candidate-peak list and the
two-group alignment map (control = denominator, target = numerator) that
findDiffPeak.DESeq2.r expects. Peak/fragment paths are written relative to the
yml's own directory, so the caller can cd into the group folder before running
the R script (which resolves inputs and writes outputs from that cwd).
"""
import argparse, os, sys
import pandas as pd
import yaml


def main():
	ap = argparse.ArgumentParser(
		prog="chip.make_group_peak_yml.py",
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="""\
Write a findDiffPeak.DESeq2.yml config for group-wise differential
peak calling.

Given a Group name and a sample.tsv, this collects every sample in
that group and builds the two-group alignment map that
findDiffPeak.DESeq2.r expects:
  - Control group = denominator, taken from the Ctrl column
                    (entries equal to 'NULL' are skipped)
  - Target group  = numerator, taken from the Name column
The candidate peak list is the per-sample HomerPeak factor peaks
(usually called without control, peak.exBL.1rpm.bed) for the target samples.

Both groups need at least 2 samples for DESeq2 dispersion
estimation; otherwise the script exits with an error.

Peak and fragment paths are written RELATIVE to the output yml's
own directory, so the run rule can cd into that folder before
invoking the R script (which resolves inputs and writes outputs
from the current working directory).""",
		epilog="""\
sample.tsv columns used:
  Name    sample id; peak/fragment files live under
          <sample-dir>/<Name>/
  Group   group label matched against --group
  Ctrl    control (input) sample id, or 'NULL' for none

example:
  chip.make_group_peak_yml.py \\
      --group H3K9me3 \\
      --sample-tsv sample.tsv \\
      --sample-dir results/samples \\
      --genome mm10 \\
      -o results/diffpeak/H3K9me3/findDiffPeak.DESeq2.yml""")
	ap.add_argument("--group", required=True,
		help="Group label to build the config for; matched against the "
			 "Group column of --sample-tsv.")
	ap.add_argument("--sample-tsv", required=True,
		help="Tab-separated sample sheet with Name, Group and Ctrl columns "
			 "('#' lines are treated as comments).")
	ap.add_argument("--sample-dir", required=True,
		help="Directory holding the per-sample output folders "
			 "(<sample-dir>/<Name>/fragment.bed.gz and HomerPeak peaks).")
	ap.add_argument("--genome", required=True,
		help="Genome build written into the config (e.g. mm10, hg38).")
	ap.add_argument("-o", "--out", required=True,
		help="Path of the yml file to write; parent directories are created "
			 "as needed. Peak/fragment paths are stored relative to this file.")

	# No arguments: show the full help instead of an argparse usage error.
	if len(sys.argv) == 1:
		ap.print_help()
		sys.exit(0)

	a = ap.parse_args()

	s = pd.read_csv(a.sample_tsv, sep="\t", comment="#", na_filter=False)
	grp = s[s.Group == a.group]
	if grp.empty:
		sys.exit("No samples for group '%s'" % a.group)

	targets = grp.Name.tolist()
	ctrls   = [c for c in dict.fromkeys(grp.Ctrl.tolist()) if c.upper() != "NULL"]
	if len(targets) < 2 or len(ctrls) < 2:
		sys.exit("Group '%s' needs >=2 target and >=2 control samples "
				 "(got %d target, %d control)" % (a.group, len(targets), len(ctrls)))

	# Paths are written relative to the yml's directory, since the run rule
	# cd's into that folder before invoking findDiffPeak.DESeq2.r
	out_dir = os.path.dirname(os.path.abspath(a.out))
	rel = lambda p: os.path.relpath(os.path.abspath(p), out_dir)
	frag = lambda n: rel(os.path.join(a.sample_dir, n, "fragment.bed.gz"))
	peak = lambda n: rel(os.path.join(
		a.sample_dir, n, "HomerPeak.factor.wo_ctrl", "peak.exBL.1rpm.bed"))

	config = {
		"pipeline": "findDiffPeak.DESeq2.r",
		"title":    a.group,
		"peak":     [peak(n) for n in targets],
		"alignment": {                              # group 1 = control (denominator)
			"Control": {n: frag(n) for n in ctrls},
			a.group:   {n: frag(n) for n in targets},   # group 2 = target (numerator)
		},
		"FC": 2, "FDR": 0.05, "PValue": "NULL",
		"independentFiltering": "FALSE",
		"ignore_ttc": "FALSE",
		"genome": a.genome,
	}
	os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
	with open(a.out, "w") as fh:
		yaml.safe_dump(config, fh, default_flow_style=False, sort_keys=False)


if __name__ == "__main__":
	main()
