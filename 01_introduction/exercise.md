# Exercises - Introduction to Python
SS 2015
HU Berlin


1. Play a round of `learn_python_light`!

    1. Start `IPython`
    2. Type `from learn_python_light import *` and follow the instructions
       (You should normally never use `from something import *`, because it messes with your namespace. But here and now: do it!)

2. Use git to track the file and changes you will produce in exercise 3.
    
    1. Create a file called `guess_a_number_yournameinshort.py` (like `guess_a_number_ka.py`). You will write the code from the next exercise into this file.
    2. Add and commit it to your local repo (from the git-exercise this morning).
    3. While writing the game from exercise 3, please commit it from time to time.
    4. Push your local commits to your remote repository.

3. Program your own game: Guess a number!

    1. Use the file you created in exercise 2.1 to write your code in. You can test the file in IPython: `run your_filename.py`
    2. Write code, that generates a random number between 0 and 100 and asks the user for a guess of that number. If the first guess isn't right, it should give a hint and ask again. And so on. To generate a random integer, you can use `randint` from the module `random`. Import that with `from random import randint`. Find out, what it does and how to use it with the help of IPython.
    3. The programm should tell the user, how many times she/he had to guess.
    4. After a round of the game is finished, the user should be asked, if she/he wants to play again.
    5. Make sure, that the user input is correct. (Has to be between 0 and 100 and should not raise an exception.)
    6. Write a 'demo version', where the computer guesses the number.

4. Learn more git!

    1. Pull changes made to the repository fugufisch/hu_bp_python_course into your fork:
       - `git remote add --track master upstream https://github.com/fugufisch/hu_bp_python_course.git`
       - `git fetch upstream`
       - `git merge upstream/master`
       - `git push origin master`

    2. If there is still time left, play with git! 
       - Add your guess_a_number to your neighbours repo (after you forked and cloned it) and issue a pull request.
       - Try branching.
       - Can you also clone your neighbours branch?
