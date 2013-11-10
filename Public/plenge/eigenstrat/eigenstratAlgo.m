Eigenstrat algorithm
Robert Chen

disease simulation:
for each SNP:
    %popAlleleFreq (p) = rand number in [0.1, 0.9]
    
    populatino allele frequency drawn from beta distribution
        parameters: p*(1-F_st) / F_st; (1-p)*(1-F_st) %mean = p, variance = F_st*p*(1-p)

RISK MODEL:
    from population l, with populatino allele frequency p_l:
    create matrix with i rows, j columns;
    i = number of SNPS;
    j = number of indivs;
    
    1 matrix for controls;
    1 matrix for cases;
         
    controls: 
        prob(genotype 0) = (1-p_l)^2;
        prob(genotype 1) = (2*p_l)*(1-p_l);
        prob(genotype 2) = (p_l)^2;

    cases:     
        prob(genotype 0) = (1-p_l)^2;
        prob(genotype 1) = (2*R*p_l)*(1-p_l);
        prob(genotype 2) = R^2 * (p_l)^2;
    
    divide all probabilities for controls/cases by (1-p_1)^2 + 2*R*p_1*(1-p_l) + R^2 * p_l^2
	
%%*RISK MODEL FOR ADMIXED POPULATION!!    
  
%%%%


INFERENCE OF AXES OF VARIATION:

    % do this for the control matrix and for case matrix:
    M = number of genotypes;
    N = number of individuals;
    
    matrix G: dimensions i (for each snp), j (for each individual);
    matrix X: modified from matrix G;
    
    for each row i in G: 
        rowSum = sum(G(i, :));
        rowMean  = mean(G(i, :)); 
     
        Define the entries of matrix X:
        for each j in row i:
            ASSIGN the following modifications of matrix G to matrix X:
                subtract rowMean from j;
            %% sum of each row should now equal 0;  
                %p_i = "posterior estimate of the unobserved underlying allele frequencyof SNP i" <<<<
                p_i = (1 + rowSum) / (2 + 2*N)
                divide j by sqrt(1 - p_i); 
        end
    end
    
    
    %%COVARIANCE MATRIX
    MATRIX psi = covariance matrix dimensions NxN;
        psi(j,j) = covariance covariance of column j of G, and column j of X;
        find the eigenvalues, eigenvectors of psi;
        eigenvals: array dimension N;
        eigenvecs: matrix dimension NxN; each ROW represents an axis of variation; %%this is essentially the ancestry matrix
        sort eigenvalues of psi;  
        %%kth axis of variation = eigenvector with the kth largest
        %%eigenvalue
        %%NOTE: k ranges from 1 to N?; imaginary eigenvals?
        
    create matrix of ancestry %%ancestry = A(j,k) = ancestry of individual j along the kth axis of variation = coordinate j of the kth eigenvector
    matrix A with dimension NxN;
       
    %ANCESTRY (jth indiv)
    A(j,k) = eigenvecs(j, k);
    
    %DECOMPOSITION  
    matrix X = U*S*transpose(V);
        %U = MxN matrix, kth col has coordinates of each SNP on kth
        %principal component
        %S = diag matrix of singular values <<<<
        %V = NxN matrix; col k has ancestries A(j,k) of each indiv j on kth
        %component
        
        
        
        
ARMITAGE TREND CHI-SQUARE STATISTIC
   armitage_chi_square = N * (squared correlation between genotype and phenotype)
         %OR (N-1) * (squared correlation) ??

GENOME-WIDE CHI-SQUARE INFLATION FACTOR FOR GENOMIC CONTROL
    median_chi_square / 0.456;  %why 0.456?
    

ADJUSTMENT OF GENOTYPES AND PHENOTYPES USING AXES OF VARIATION
    recall: G(i,j) = genotype of individual j;
            A(j,k) = ancestry of individual along the kth axis of var %one ancestry value per SNP
    %gamma = regression coefficient; one value for each SNP
    gammma = matrix with dim NxN;
    for a from 1 to N: %each row
        for b from 1 to N: %each col value in that row (each person)
            ancestry = ???;
            gammaNumerVal = ancestry*g(a,b);
            gammaDenomVal = ancestry^2;
            gamma = gammaNumerVal / gammaDenomVal;
        end
    end
    
    for each Gadjusted: = Gadjusted(i,j) = G(i,j) - gamma(i) * a(j);
    %%repeat same process for adjusted phenotype matrix
    
    
        
    Gadjusted (i,j) = G(i,j) - gamma(i)*A(j);
        
            
    
        
        
        
    
    
                
                
    
            
            
        
    
