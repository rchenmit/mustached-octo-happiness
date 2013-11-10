Eigenstrat algorithm

disease simulation:
for each SNP:
    %popAlleleFreq (p) = rand number in [0.1, 0.9]
    populatino allele frequency drawn from random distribution
        parameters: p*(1-F_st) / F_st; (1-p)*(1-F_st) %mean = p, variance = F_st*p*(1-p)

RISK MODEL:
    from population l, with populatino allele frequency p_l:
    
    controls: 
        prob(genotype 0) = (1-p_l)^2;
        prob(genotype 1) = (2*p_l)*(1-p_l);
        prob(genotype 2) = (p_l)^2;

    cases:     
        
