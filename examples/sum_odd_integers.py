import logging

# Use logging rather then print statements to track workflows and record exceptions
logging.basicConfig(level=logging.INFO,  # DEBUG, INFO, WARN, ERROR
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.FileHandler('sum_odd_integers.log'),
                              logging.StreamHandler()],)


"""
Given two positive integers a and b, return the sum of all odd integers from a to b inclusively

So if a = 50 and b = 100 the output = 1875
Challenge: what is the output if a = 10 and b = 25?
"""


def sum_odd_integers(**kwargs):
    """
    Function that sums the odd integers in a defined ranbe
    :param qa:
    :param qb:
    :return:
    """
    logging.info("sum odd integers between {} and {}".format(str(kwargs["a"]), str(kwargs["b"] + 1)))
    sum_up = 0
    for i in range(kwargs["a"], kwargs["b"] + 1):
        logging.debug("Is {} an odd number? If not pass".format(str(i)))
        if i % 2 != 0:  # Modulo operator.
            logging.info("add = {} to {}".format(str(i), str(sum_up)))
            sum_up += i
    return sum_up


if __name__ == '__main__':
    """
    This code allows me to directly call the script and use a test string e.g. the one below
    """
    odd_integers_summed = sum_odd_integers(a=50, b=100)
    print(odd_integers_summed)
