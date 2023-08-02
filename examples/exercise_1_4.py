import logging

# Use logging rather than print statements to track workflows and record exceptions
logging.basicConfig(level=logging.DEBUG,  # DEBUG, INFO, WARN, ERROR
                    format="%(asctime)s [%(levelname)s] %(message)s",
                    handlers=[logging.FileHandler('chunking.log'),
                              logging.StreamHandler()],)


def chunk_string(query_sequence, chuk_by):
    """
    Function to split sequence string into chunks of a value
    :param query_sequence: DNA sequence
    :param chuk_by: Integer. Chunk by
    :return:
    """
    logging.info("Chunk {} into blocks of {}".format(query_sequence, str(chuk_by)))

    my_list = []
    while query_sequence:
        my_list.append(query_sequence[:chuk_by])
        query_sequence = query_sequence[chuk_by:]
    return " ".join(my_list)


if __name__ == "__main__":
    string = "aggagtaagcccttgcaactggaaatacacccattg"
    chunk_length = 5
    print(chunk_string(string, chunk_length))
