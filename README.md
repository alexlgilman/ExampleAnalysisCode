# ExampleAnalysisCode
Example BESIII Analysis code
- DStarIDAlg: BESIII Analysis code package that identifies D{\*+} candidates through the decay chain D{*+}->Pi+ D^{0}, D^{0}->K-Pi.
- DStarID_ana.txt: BOSS Job Options file to run the compiled code.
- MCGeneration: Shows example of how to generate MC samples within BOSS

#How to Write a BOSS script
First, setup your BOSS environment according to https://docbes3.ihep.ac.cn/~offlinesoftware/index.php/How_to_setup_BOSS_environment_on_CentOS7. To setup BOSS after the first setup, you’ll just need to run the setup.sh files in your cmthome directory and the TestRelease directory.
To create an analysis script, go to your workarea directory and run “cmt create ScriptName VersionNumber” where VersionNumber is in the format XX-XX-XX
“cd ScriptName/VersionNumber”
“mkdir share”  and place a share file in this directory. An example sharefile is included in the attachment
“mkdir ScriptName;” for the header file. Example is included in the attachment.
Modify the file “cmt/requirements” to include the appropriate BOSS packages.
Edit the source files in the “src” directory. Some directories will need to be updated.
To compile, cd to the “cmt” directory and “cmt br make”
To run the code, first source “ScriptName/VersionNumber/cmt/setup.sh” and then run “boss.exe JobOptionsFile.txt”. An example JobOptionsFile.txt titled DStarID_ana.txt is included in the attachement. Information about the MC samples is at https://docbes3.ihep.ac.cn/~charmgroup/index.php/MC_Samples.
