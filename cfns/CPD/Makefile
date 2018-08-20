INSTALLPATH = /home/tschiex/Dev/msd/positive

PYTHON = python3
PYPATH = $(INSTALLPATH)/exes
PARPATH = ./pars
VPATH = $(PYPATH):$(PARPATH)

ALLPDBS = $(wildcard *.pdb)
ALLOPTS = $(filter-out positive.opt, $(ALLPDBS:.pdb=.opt))
ALLENUMS = $(filter-out positive.enum, $(ALLPDBS:.pdb=.enum))
ALLGMECS =  $(filter-out positive.gmec, $(ALLPDBS:.pdb=.gmec))
ALLCFNS =  $(filter-out positive.cfn.gz, $(ALLPDBS:.pdb=.cfn.gz))

NJOBS = $(shell python -c 'import multiprocessing as mp; print(mp.cpu_count()+1)')

NO_COLOR=\033[0m
OK_COLOR=\033[32;01m

all : $(ALLOPTS) $(ALLENUMS)

#$(info $$ALLCFNS is [${ALLCFNS}])

.PHONY: clean superclean
clean:
	\rm -f *.rlx *.cfn.gz *.conf *.gmec *.pwa *.show *.enum *.tb2 *.opt
	\rm -f *.gmec+ *.enum+ *.opt+ *.pwa+ *.show+ *.wlst

superclean: clean
	\rm -rf $(PYPATH)/tb2 $(PYPATH)/toulbar2

tb2:
	@printf "$(OK_COLOR)Cloning toulbar2 $(NO_COLOR)\n"
	git clone http://mulcyber.toulouse.inra.fr/anonscm/git/toulbar2/toulbar2.git $(PYPATH)/tb2

toulbar2: tb2 -lgmp -lboost_graph -lboost_iostreams 
	@printf "$(OK_COLOR)Compiling toulbar2 $(NO_COLOR)\n"
	cd $(PYPATH)/tb2/toulbar2 && \
	git checkout cpd && \
	git checkout `git describe --abbrev=0` && \
	mkdir -p build && \
	cd build && \
	cmake .. && \
	make -j$(NJOBS)
	cd $(PYPATH) && \
	ln -s ./tb2/toulbar2/build/bin/*/toulbar2 toulbar2

showpars:
	@printf "Fixed seed: "
	@cat pars/fixedseed
	@printf "Native bias: "
	@cat pars/nativebias
	@printf "PSSM bias: "
	@cat pars/pssmbias
	@printf "Number of fixbb: "
	@cat pars/nbfixbb
	@printf "Number of fastrelax: "
	@cat pars/nbrelax
	@printf "Precision: "
	@cat pars/precision
	@printf "Rotamers: "
	@cat pars/rotamers
	@printf "Score function: "
	@cat pars/scorefunction
	@printf "Energy visualization: "
	@cat pars/energy_vis
	@printf "H-bonds visualization: "
	@cat pars/hbonds_vis
	@printf "Enumeration threshold: "
	@cat pars/threshold
	@printf "Verbosity: "
	@cat pars/verbose

pymolserver:
	@touch  ~/.pymolrc
	@if grep PyMOL-RosettaServer.py ~/.pymolrc > /dev/null 2>&1 ; then \
	echo "Already installed";\
	else \
	echo "run $(PYPATH)/PyMOL-RosettaServer.py" >>  ~/.pymolrc ;\
	echo "~/.pymolrc startup script modified"; \
	fi

positive.cfn.gz positive.wlst: $(ALLCFNS) $(PYPATH)/make-pos.py
	@printf "$(OK_COLOR)Generating positive multistate matrix $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/make-pos.py $(ALLCFNS) 

positive.gmec positive.tb2: positive.cfn.gz toulbar2 $(PYPATH)/tb2cpd.py
	@printf "$(OK_COLOR)Looking for optimal positive multistate sequence $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/tb2cpd.py --dogmeconf -i $<

positive.enum: positive.cfn.gz positive.tb2 toulbar2 threshold sampling $(PYPATH)/tb2cpd.py
	@printf "$(OK_COLOR)Enumerating suboptimal positive multistate sequences $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/tb2cpd.py --doenum -S -i $<

%.gmec+: positive.gmec positive.wlst $(PYPATH)/proj-pos.py
	@printf "$(OK_COLOR)Projecting multistate optimum sequence/conformation $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/proj-pos.py --conf  $<

%.enum+: positive.enum positive.wlst $(PYPATH)/proj-pos.py
	@printf "$(OK_COLOR)Projecting multistate suboptimal sequences $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/proj-pos.py --enum  $<	

%.opt+ %.show+: %.gmec+ %.rlx %.resfile scorefunction rotamers $(PYPATH)/tb2cpd.py
	@printf "$(OK_COLOR)Building the optimum multistate pose for $< $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/tb2cpd.py --dogmepose --positive -i $<

%.pwa+: %.gmec+ %.rlx %.resfile scorefunction rotamers $(PYPATH)/tb2cpd.py
	@printf "$(OK_COLOR)Computing multistate pairwise energies and distances for $< $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/tb2cpd.py --dopairwise --positive -i $<

pars:
	cp -rf $(INSTALLPATH)/$(PARPATH) .

scorefunction nbrelax precision rotamers psimatrix nativebias nbfixbb: pars

%.rlx: %.pdb scorefunction nbrelax $(PYPATH)/tb2cpd.py
	@printf "$(OK_COLOR)Relaxing $< $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/tb2cpd.py --dorelax -i $< 

%.resfile:
	@printf "$(OK_COLOR)No resfile found, copying default resfile to $@ $(NO_COLOR)\n"
	cp $(PARPATH)/resfile $@

%.cfn.gz: %.rlx %.resfile scorefunction precision rotamers psimatrix nativebias $(PYPATH)/tb2cpd.py
	@printf "$(OK_COLOR)Computing pairwise energy design matrix for $< $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/tb2cpd.py --doematrix -i $< 

%.conf: %.rlx %.resfile scorefunction rotamers nbfixbb $(PYPATH)/tb2cpd.py
	@printf "$(OK_COLOR)Computing upper bounding conformation for $< $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/tb2cpd.py --dofixbb -i $< 

%.tb2 %.gmec: %.cfn.gz toulbar2 $(PYPATH)/tb2cpd.py | %.conf
	@printf "$(OK_COLOR)Computing optimal sequence/conformation for $< $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/tb2cpd.py --dogmeconf -i $<

%.opt %.show: %.gmec %.rlx %.resfile scorefunction rotamers $(PYPATH)/tb2cpd.py
	@printf "$(OK_COLOR)Building the optimum pose for $< $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/tb2cpd.py --dogmepose -i $<

%.pwa: %.gmec %.rlx %.resfile scorefunction rotamers $(PYPATH)/tb2cpd.py
	@printf "$(OK_COLOR)Computing pairwise energies and distance for $< $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/tb2cpd.py --dopairwise -i $<

%.enum: %.cfn.gz %.tb2 toulbar2 threshold sampling $(PYPATH)/tb2cpd.py
	@printf "$(OK_COLOR)Enumerating suboptimal sequences for $< $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/tb2cpd.py --doenum -S -i $<

%.scan: %.cfn.gz %.tb2 toulbar2 $(PYPATH)/tb2cpd.py
	@printf "$(OK_COLOR)Scanning for less mutations for $< $(NO_COLOR)\n"
	$(PYTHON) $(PYPATH)/tb2cpd.py --doscan -S -i $<



.PRECIOUS: %.rlx %.cfn.gz %.scan %.conf %.gmec %.enum *.pwa %.tb2 %.opt %.show %.resfile

