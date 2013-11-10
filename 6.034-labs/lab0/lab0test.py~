#lab 0 test


def cube(x):
    return x*x*x

def factorial(x):
    return reduce(lambda i,j: i*j, range(1,x+1))

def count_pattern(pattern, lst):
	currentSearch = 0
	appearCount = 0
	for i in range(len(lst)):
		gap = 0	
		currentSearchIndex = i
		for j in range(len(pattern)):
			if lst[currentSearchIndex] == pattern[j]:
				currentSearchIndex += 1
				gap += 1 
				if j == (len(pattern)-1):
					appearCount += 1
			else:
				break
	return appearCount	
