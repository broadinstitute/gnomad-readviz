from typing import List


def get_sample_ids(ids_file: str, header: bool = False) -> List[str]:
    """
    Open file with sample IDs and stores IDs in a list.

    :param str ids_file: Path to file containing sample IDs.
    :param bool header: Whether IDs file has a header line. Default is False.
    :return: List of sample IDs
    :rtype: List[str]
    """
    sample_ids = []
    with open(ids_file) as i:
        if header:
            header = i.readline()
        for line in i:
            sample_ids.append(line.strip())
    return sample_ids
