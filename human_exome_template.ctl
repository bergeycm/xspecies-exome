GENERAL-INFO-START

	seq-file            SEQ_FILE_HERE
	trace-file          SEQ_FILE_HERE.trace.log
	locus-mut-rate      VAR 1.0

	burn-in             0
	mcmc-iterations	    5000000
	iterations-per-log  1000
	logs-per-line       100
	mcmc-sample-skip   10
	
#	tau-theta-alpha		1.0
#	tau-theta-beta		10000.0
	tau-theta-print     10000.0

#	find-finetunes		TRUE
	finetune-coal-time	0.3		
	finetune-mig-time	0.3		
	finetune-mig-rate	0.02
	
	finetune-theta      0.04
	finetune-tau        0.0000008

	finetune-locus-rate	0.5
	finetune-mixing	0.003

##	mig-rate-alpha		0.002
##	mig-rate-beta		0.0000001

GENERAL-INFO-END

CURRENT-POPS-START	

	POP-START
		name		EUROPEAN
		samples		venter d na12891 d
		theta-alpha	1.0
		theta-beta      10000
	POP-END

	POP-START
		name		AFRICAN_NOT_BUSHMAN
		samples		na18507 d
		theta-alpha	1.0
		theta-beta      10000
	POP-END
	
	POP-START
		name		CHINESE
		samples		hanChinese d
		theta-alpha	1.0
		theta-beta      10000
	POP-END

	POP-START
		name		KOREAN
		samples		sjk d
		theta-alpha	1.0
		theta-beta      10000
	POP-END

	POP-START
		name		BUSHMAN
		samples		kb1 d
		theta-alpha	1.0
		theta-beta      10000
	POP-END

	POP-START
		name		CHIMP
		samples		chimp h
		theta-alpha	1.0
		theta-beta      10000
	POP-END
	
CURRENT-POPS-END

ANCESTRAL-POPS-START

	POP-START
		name			ASIAN
		children		CHINESE KOREAN
		tau-alpha	1.0
		tau-beta	30000
		theta-alpha	1.0
		theta-beta	10000
		tau-initial	0.0000333
	POP-END

	POP-START
		name			EURASIAN
		children		EUROPEAN ASIAN
		tau-alpha	1.0
		tau-beta	30000
		theta-alpha	1.0
		theta-beta	10000
		tau-initial	0.0000333
	POP-END	

	POP-START
		name			NOT_BUSHMAN
		children		EURASIAN AFRICAN_NOT_BUSHMAN
		tau-alpha	1.0
		tau-beta	25000
		theta-alpha	1.0
		theta-beta	10000
		tau-initial	0.00004
	POP-END	

	POP-START
		name			HUMAN
		children		NOT_BUSHMAN BUSHMAN
		tau-alpha	1.0
		tau-beta	10000
		theta-alpha	1.0
		theta-beta	10000
		tau-initial	0.0001
	POP-END	

	POP-START
		name			ROOT
		children		HUMAN CHIMP
		tau-alpha	1.0
		tau-beta	1000
		theta-alpha	1.0
		theta-beta	10000
		tau-initial	0.001
	POP-END

ANCESTRAL-POPS-END
