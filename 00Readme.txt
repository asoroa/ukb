
To create a branch

1) make sure all changes are commited

2) create a directory under Mcr_main (via emacs), or

% svn mkdir svn+ssh://siuc04/sc02a1/users/ccpsoeta/svnRep/ukb/Main_tags/s2aw_exp

3) make a copy (branch) of current HEAD directories Main, Preproc to new created directory

% svn copy svn+ssh://siuc04/sc02a1/users/ccpsoeta/svnRep/ukb/Main svn+ssh://siuc04/sc02a1/users/ccpsoeta/svnRep/ukb/Main_tags/s2aw_exp
% svn copy svn+ssh://siuc04/sc02a1/users/ccpsoeta/svnRep/ukb/Preproc svn+ssh://siuc04/sc02a1/users/ccpsoeta/svnRep/ukb/Main_tags/s2aw_exp


///////////////////////////
graphviz

neato w.dot -Tsvg > w.svg
gqview w.svg 

