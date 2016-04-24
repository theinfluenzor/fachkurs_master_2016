import string
import random

passwords = {}
pw_length = 8
allowed_characters = string.ascii_letters + string.punctuation  + string.digits

def generate_pw(n, allowed):
    pw = ''
    for i in range(pw_length):
        pw += random.choice(allowed)
    return pw

selection = input("(r)etrieve or (g)enerate password? ")
if selection == 'r':
    print(passwords.keys())
    print(passwords[input("Which password would you like? ")])
elif selection == 'g':
    passwords[input("Enter account name: ")] = generate_pw(pw_length, allowed_characters)
