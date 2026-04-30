import math


def GetData(filename):
    frequencies = {}
    with open(filename) as f:
        data = f.read()
    lines = [s.strip().split() for s in data.splitlines()]
    lines.pop(0)
    lines.pop(0)
    for line in lines:
            frequencies[line[0]]= line[3]
    
    return frequencies


def DotProduct(vector1, vector2):
    dot_product = 0
    for key in vector2:
        if key in vector1:
            dot_product +=  (float(vector1[key]) * float(vector2[key]))
    return dot_product



def VectorMagnitude(vector):
    magnitude = 0
    for key in vector:
        magnitude = magnitude + (float(vector[key]) ** 2)
    magnitude = math.sqrt(magnitude)
    return magnitude


def FindAngle(dot_product, magnitude1, magnitude2):
    result = dot_product / (magnitude1 * magnitude2)
    return result


vector1 = GetData('WT_singlebase_deletionsv2_Trinuc_trinucmutation_strand_background.txt')
vector1b = GetData('rad30_singlebase_deletions_Trinuc_trinucmutation_strand_background.txt')

vector2 = GetData('WT_singlebase_deletions_Trinuc_trinucmutation_strand_background.txt')
vector2b = GetData('rad16_singlebase_deletions_Trinuc_trinucmutation_strand_background.txt')

vector3 = GetData('WT_singlebase_deletions_Trinuc_trinucmutation_strand_background.txt')
vector3b = GetData('rad26_singlebase_deletions_Trinuc_trinucmutation_strand_background.txt')


dot_product = DotProduct(vector1, vector1b)
magnitude1 = VectorMagnitude(vector1)
magnitude2 = VectorMagnitude(vector1b)
cosine = FindAngle(dot_product, magnitude1, magnitude2)
print(cosine)

dot_productb = DotProduct(vector2, vector2b)
magnitude1b = VectorMagnitude(vector2)
magnitude2b = VectorMagnitude(vector2b)
cosineb = FindAngle(dot_productb, magnitude1b, magnitude2b)
print(cosineb)

dot_productc = DotProduct(vector3, vector3b)
magnitude1c = VectorMagnitude(vector3)
magnitude2c = VectorMagnitude(vector3b)
cosinec = FindAngle(dot_productc, magnitude1c, magnitude2c)
print(cosinec)
