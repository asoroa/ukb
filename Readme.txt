UKB: Graph Based Word Sense Disambiguation and Similarity

UKB is a collection of programs for performing graph-based Word Sense
Disambiguation and lexical similarity/relatedness using a pre-existing
knowledge base.

UKB has been developed by the IXA group in the University of the Basque
Country. UKB applies the so-called Personalized PageRank on a Lexical
Knowledge Base (LKB) to rank the vertices of the LKB and thus perform
disambiguation. The details of the method are described in [1]. It has also
been applied on WSD on specific domains [2,5] and Named Entity
Disambiguation [6]. The algorithm can also be used to calculate lexical
similarity/relatedness of words/sentences. See [3,4,6] for applications of
UKB to similarity.

Visit http://ixa2.si.ehu.es/ukb/ for more information about UKB.

The latest source code for git can be found here:

https://github.com/asoroa/ukb.git


References
**********

[1] Eneko Agirre and Aitor Soroa. 2009. Personalizing PageRank for
Word Sense Disambiguation. Proceedings of the 12th conference of the
European chapter of the Association for Computational Linguistics
(EACL-2009). Athens, Greece.

[2] Eneko Agirre, Oier Lopez de Lacalle and Aitor
Soroa. 2009. Knowledge-based WSD and specific domains: performing over
supervised WSD. Proceedings of IJCAI. Pasadena, USA.

[3] Eneko Agirre, Enrique Alfonseca, Keith Hall, Jana Kravalova,
Marius Pasca and Aitor Soroa. 2009. A Study on Similarity and
Relatedness Using Distributional and WordNet-based
Approaches. Proceedings of NAACL-HLT 09. Boulder, USA.

[4] Eneko Agirre, Montse Cuadros, German Rigau and Aitor Soroa. 2010.
Exploring Knowledge Bases for Similarity. Proceedings of LREC
2010. Valletta, Malta.

[5] Eneko Agirre, Aitor Soroa, Mark Stevenson. 2010. Graph-based Word
Sense Disambiguation of Biomedical Documents. Bioinformatics, Oxford
University Press. Bioinformatics Vol. 26(22) pp: 2889-2896

[6] Eneko Agirre, Ander Barrena and Aitor Soroa. 2015. Studying the
Wikipedia Hyperlink Graph for Relatedness and Disambiguation. 
http://arxiv.org/abs/1503.01655


Files under this catalogue:
***************************

src/                    Source code. See README, INSTALL, LICENSE.
bin/                    Statically compiled binaries for x86-64 linux
                        platforms.
UKBsim/                 Scripts for similarity. See README, INSTALL, LICENSE.
scripts/                Scripts for converting Wordnet to ukb input files

Check README files in the respective catalogue.

Getting the resources
*********************

* WordNet

- English WordNet (versions 1.7 and 3.0) can be found here:
http://ixa2.si.ehu.es/ukb/lkb_sources.tar.bz2

- Spanish WordNet is here:
http://ixa2.si.ehu.es/ukb/lkb_sources.tar.bz2

- More WordNet versions can be converted to UKB format using the scripts in
  "scripts/" directory.

* Wikipedia

- A graph derived from Wikipedia (dump 04 April 2013) is here:
http://ixa2.si.ehu.eus/ukb/ukb-wiki.tar.bz2
