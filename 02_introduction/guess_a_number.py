from random import randint


def guess_a_number():
    """Game to guess a number the computer randomly generated."""

    random_number = randint(0, 100)

    demo_or_play = unicode(raw_input("Do you want the computer to guess the number? If so, say 'demo'! "))
    if demo_or_play == 'demo':
        demo_a_number(random_number)
    else:
        guess = check_raw("Please enter an integer between 0 and 100: ")
        count_guesses = 1
        while guess != random_number:
            count_guesses = count_guesses + 1
            guess = evaluate_my_number(guess, random_number)
        else:
            print('Hooray! Your guess is right!\nYou needed {} guesses.'.format(count_guesses))

    new_game = unicode(raw_input("Do you want to play again? If so, say 'yes'! "))
    if new_game == 'yes':
        guess_a_number()
    else:
        print('Goodbye!')


def evaluate_my_number(guess, random_number):
    """Is the guess to high or to low? Guess again!"""

    if guess < random_number:
        print 'Too low!',
    else:
        print 'Too high!',
   
    guess = check_raw()
    return guess


def check_raw(print_string='Please try again: '):
    """Gets the string, raw_input should print, checks and returns the input."""

    try:
        checked_int = int(raw_input(print_string))
        if checked_int < 0 or checked_int > 100:
            print 'Your number has to be between 0 and 100!',
            checked_int = check_raw()
    except ValueError:
        print 'That was not an integer!',
        checked_int = check_raw()

    return checked_int


def demo_a_number(random_number):
    """The computer tries to guess the number"""

    current_number = 50
    lower_bound = 0
    upper_bound = 100

    count_computer_guesses = 1
    while current_number != random_number:
        count_computer_guesses = count_computer_guesses + 1
        print('The computer guessed {}'.format(current_number))
        if current_number < random_number:
            print('That was too low.')
            lower_bound = current_number
        else:
            print('That was too low.')
            upper_bound = current_number
        current_number = (lower_bound + upper_bound) / 2
    else:
        print('The computer guessed {}\nThat was right! \
        After {} guesses'.format(current_number, count_computer_guesses))


guess_a_number()
