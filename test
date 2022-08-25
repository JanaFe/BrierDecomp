#### Testing BrierDecomp Function with two tests:

# TEST 1: with number of bins specified
p = np.array([0.1923, 0.9532, 0.8834, 0.2420, 0.0521, 0.6787, 0.4321]) 
y = np.array([0,1,1,0,0,0,1])
bins = 10
bias_corrected=False

est, sd = BrierDecomp(p,y, bins=10, bias_corrected=False)
print(est)
print(sd)



# TEST 2: With bins specified and bias correction
p = np.array([0.94725436,0.9053663,0.15559863,0.18543585,0.24392618,0.97578883,0.13123178,0.09639295,0.16112317,0.115615964,
              0.7339622,0.22357145,0.33462253,0.17155208,0.14778529,0.09295429,0.47029933,0.7423472,0.3012784,0.17553441]) 
y = np.array([1,1,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,0])
bins = np.quantile(p, q=np.arange(0.0, 1.1, 0.1)) #until 1.1, so that 1.0 is still included
bias_corrected=False

est, sd = BrierDecomp(p,y, bins=bins, bias_corrected=True)
print(est)
print(sd)


#Results matched the output of BrierDecomp in R. 
