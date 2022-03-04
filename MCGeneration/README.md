*--------------*
Example MC Generation Scripts using BesEvtGen
Author: Alex Gilman
Last updated 25-05-21
*--------------*

Included in this directory are two .DEC or "decay" files which prescribe the allowed decays for included particles. All other particle decays are inherited from the default .DEC file in your BOSS version.

For more info on BesEvtGen see https://docbes3.ihep.ac.cn/cgi-bin/DocDB/ShowDocument?docid=102
For available decay models see https://docbes3.ihep.ac.cn/cgi-bin/DocDB/ShowDocument?docid=162

The file sim.txt simulates simulates the fundamental physics and outputs a .rtraw file. A decay file must be supplied, as well as a random number to seed MC generation.

The file rec.tx takes a .rtraw file as input and outputs a .dst file with the reconstructed information.

