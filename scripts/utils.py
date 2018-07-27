def exclude_sample(sample, blacklist):
    for pattern in blacklist:
        if pattern in sample:
            return True
    return False
