#!/usr/bin/env python3
#this script prints the numbers between 0 and 99 no divisible by 7
nums = list(range(100))
for n in nums:
    if n == 0 or n%7 > 0:
        print(n)