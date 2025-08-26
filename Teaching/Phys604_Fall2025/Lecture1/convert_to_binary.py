
# Input decimal and number of digits of precision
decimal = float(input("Enter a float to convert: "))
prec_digits=int(input("Enter an integer for precision: "))

# Convert to binary
int_version=int(decimal*2**prec_digits)
binary = bin(int_version)
decimal_length=prec_digits-len(binary.split('b')[1])

# Output formatted binary result
print(float(binary.replace('0b',''))*(1.0*10**(-1*prec_digits)))
