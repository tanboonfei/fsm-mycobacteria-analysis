RELAX outputs reported in the manuscript.

Primary panel: 52 single-copy orthologues, 37 taxa, 63,840 nt (21,280 codons)
  alignment            three_way.codon.fa
  branch-labelled tree common_reference_FSM_tips_stem.nwk, common_reference_LEP_tips_stem.nwk
                       common_reference_FSM_full_crown.nwk, common_reference_LEP_full_crown.nwk
                       Each carries explicit {Test} and {Reference} labels; 39 reference
                       branches are labelled and HyPhy fits 38, dropping one of zero length.
  output               common_reference_FSM_tips_stem.json    K = 1.0438  Table S11 row 1
                       common_reference_FSM_full_crown.json   K = 1.0462  Table S11 row 2
                       common_reference_LEP_tips_stem.json    K = 0.7579  Table S11 row 3
                       common_reference_LEP_full_crown.json   K = 0.7709  Table S11 row 4
  The matching *.command.json files record the invocation.
  sha256 of each output is given in Table S11.

Leprosy-free panel: 489 single-copy orthologues, 34 taxa, 170,788 codons
  alignment            super488.fas
  branch-labelled tree rx489_FSM.nwk   (27 test branches, 36 fitted reference)
  output               relax488.json    K = 1.0724  Table S11g
  md5 of the output is given in Table S11g.

In every run the reference set is RGM plus SGM only. The superseded runs from the
earlier design, in which the other reduced lineage was part of the reference set,
are recorded in Tables S11e and S11f and are not used anywhere in the main text. The
trees from that earlier design (rx_*.nwk, crown_*.nwk) carry {Test} labels only, so
every unlabelled branch including the other reduced lineage would enter the reference;
they are kept in ../_superseded/relax_trees_old/ and must not be used to reproduce
the runs above.

Stem-branch test: the same 52-orthologue alignment with only the FSM stem branch
in the test set, every FSM crown and terminal branch unclassified.
  branch-labelled tree stem_only_FSM.nwk   (1 test branch, 39 labelled reference,
                                            38 fitted, 31 unclassified)
  output               relax_stem_only.json   K = 0.940, LRT 2.07, p = 0.15
  sha256 of the output is given in Table S11h.

Joint contrast and synonymous-rate-variation runs (new_runs/)
  The same 52-orthologue alignment, relabelled so that one reduced lineage is the
  test set and the other the reference.
  output               new_runs/joint_FSM_vs_LEP.json   Table S11j
                       new_runs/joint_LEP_vs_FSM.json   Table S11j
  The same four primary runs repeated under --srv Yes.
  output               new_runs/srv_FSM_tips_stem.json  Table S11k
                       new_runs/srv_LEP_tips_stem.json  Table S11k
  new_runs/run.sh records the invocations and the matching .log files the console
  output of each run.

The "file name" field in relax488.json records the alignment at its location in
this repository; that run was executed from a scratch directory. No result field
was altered, and the md5 in Table S11g is of the file as distributed.

Note: "Table S9" above is the original numbering; in the submitted Additional file 2 these RELAX tables are Table S8 (see scripts/TABLE_NUMBERING.md).
