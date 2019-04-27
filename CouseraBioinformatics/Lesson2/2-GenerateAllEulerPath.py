from collections import defaultdict
import copy

pairs = '(AG|AG), (AG|TG), (CA|CT), (CT|CA), (CT|CT), (GC|GC), (GC|GC), (GC|GC), (TG|TG)'
k = 2
d = 1
pairs = list(map(lambda x: [x[1:3], x[4:6]], pairs.split(', ')))

debruijn = defaultdict(list)
for pair in pairs:
	debruijn[pair[0][:-1] + pair[1][:-1]].append(pair[0][1:] + pair[1][1:])

cnt = defaultdict(int)
for pos1 in list(debruijn.keys()):
	cnt[pos1] += 1
	for pos2 in debruijn[pos1]:
		cnt[pos2] -= 1
for pos in cnt:
	if cnt[pos] == 1:
		start = pos
	elif cnt[pos] == -1:
		end = pos
#print(cnt)
del cnt, pos, pos1, pos2, pair
debruijn[end].append(start)
#fill the cycle

ret = []
keylist = []
for key in debruijn.keys():
	if len(debruijn[key]) > 1:
		keylist.append(key)
	elif len(debruijn[debruijn[key][0]]) > 1:
		keylist.append(key)
#print(keylist)
for node in keylist:
	temp = copy.deepcopy(debruijn)
	string = []
	nextnode = node
	stack = [node]
	while stack != []:
		if temp[nextnode]:
			nextnode = temp[nextnode].pop()
			stack.append(nextnode)
		else:
			del temp[nextnode]
			string.append(stack.pop())
			try:
				nextnode = stack[-1]
			except:
				break
	string = string[::-1][:-1]
	string = string[string.index(start):] + string[:string.index(start)]
	if string not in ret:
		ret.append(string)
print(ret)
#find all Euler Cycle

i = 0
while i < len(ret):
	string = ret[i]
	valid = 1
	for j in range(k + d, len(string)):
		if string[j][0] != string[j - k - d][1]:
			valid = 0
			break
	if valid == 0:
		del ret[i]
	else:
		ret[i] = ''
		for j in range(k + d):
			ret[i] += string[j][0]
		for j in range(len(string)):
			ret[i] += string[j][-1]
		i += 1
#find all valid cycle

print(ret)