# LocusTagMapper

Maps a list of locus tags of different species to your species of choice

created by Andreas U. Mueller, 2020


What should the program be capable of?

Batch conversion of locus tags from multiple species to one target species

What is it NOT for?

Although LTMapper can do it, it's not intended to look up just single locus tags.

This can be done easily using existing online tools & databases.


### TODO
- user-provided target/query genomes / existing file import
- generate template config file (with current settings) on the fly if it does not exist (for future runs)

### Changelog
v0.4 - 20201003
- timestamped log (to match result files)
- revised final statistics output
- fixed: lt_getid() accepted partial matches leading to falsely assigned sequences/IDs
v0.3 - 20201001
- correct interpretation of no hits BLAST results
- fixed query CLI argument
- various minor code improvements
v0.2 - 20200926
- added CLI parameters
- reconciles CLI with config file inputs (CLI takes precedence over config)
- reports more statistics
- timestamped results
v0.1 - 20200506
- initial version, basic functionality
