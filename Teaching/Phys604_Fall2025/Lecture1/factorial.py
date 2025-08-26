
in_int = input("Enter an integer: ")

out_fact=1
for loop_iter in range(1,int(in_int)+1):
    out_fact*=loop_iter

print("The factorial of ",in_int, " is ",out_fact)
print("Memory needed: ",out_fact.bit_length())
