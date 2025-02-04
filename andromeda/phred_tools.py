"""
phred_tools.py
Marcus Viscardi,    February 04, 2025

Some helper functions to deal with Phred scores and quality values.

Major use for this will be averaging the Phred scores from multiple reads to get a consensus Phred score.

Based on some stuff from the following sources:
https://www.drive5.com/usearch/manual/quality_score.html
https://stackoverflow.com/questions/58886389/how-can-i-create-a-dictionary-that-contains-all-the-phred-scores-0-93-encoded-in
https://nanoporetech.com/support/software/data-analysis/where-can-i-find-out-more-about-quality-scores
http://samtools.github.io/hts-specs/
"""
import random
from typing import List
from math import log10
from tqdm.auto import tqdm

from utils import clamp


class NucleotideQuality:
    """
    Class to handle individual Phred quality scores for nucleotides.
    """
    def __init__(self, char: str = None, q_score: int = None, prob_error: float = None):
        self.char = char
        self.q_score = q_score
        self.prob_error = prob_error
        if not any((
                (char is not None),
                (q_score is not None),
                (prob_error is not None),
        )):
            raise ValueError("At least one of the following must be set: char, q_score, prob_error.")

    def __add__(self, other):
        if isinstance(other, NucleotideQuality):
            return NucleotideQuality(prob_error=(self.to_prob_error() + other.to_prob_error()))
        elif isinstance(other, (int, float)):
            return NucleotideQuality(prob_error=(self.to_prob_error() + other))
        else:
            raise TypeError(f"Unsupported type for addition with NucleotideQuality class: {type(other)}")
    
    def __mul__(self, other):
        if isinstance(other, NucleotideQuality):
            raise NotImplementedError("Multiplication of two NucleotideQuality objects is not supported "
                                      "(and doesn't make much sense).")
        elif isinstance(other, (int, float)):
            return NucleotideQuality(prob_error=(self.to_prob_error() * other))
        else:
            raise TypeError(f"Unsupported type for multiplication with NucleotideQuality class: {type(other)}")

    def to_phred_char(self) -> str:
        if self.char is not None:
            return self.char  # Easy answer
        elif self.q_score is not None:
            return self._phred_q_score_int_to_char()  # One-step conversion
        elif self.prob_error is not None:
            self.to_q_score()  # Calculate Q-score first
            return self._phred_q_score_int_to_char()  # Then convert to character
        else:
            raise ValueError("No Phred score information set. This should be unreachable?!")

    def to_q_score(self) -> int:
        if self.q_score is not None:
            return self.q_score  # Easy answer
        elif self.char is not None:
            return self._phred_char_to_q_score_int()  # Convert character to Q-score. One-step conversion
        elif self.prob_error is not None:
            return self._prob_error_to_phred_q_score_int()  # Also one step, since Q-score is in the middle
        else:
            raise ValueError("No Phred score information set.")

    def to_prob_error(self) -> float:
        if self.prob_error is not None:
            return self.prob_error  # Easy answer
        elif self.q_score is not None:
            return self._phred_q_score_int_to_prob_error()  # Convert Q-score to probability of error
        elif self.char is not None:
            self._phred_char_to_q_score_int()  # Convert character to Q-score first
            return self._phred_q_score_int_to_prob_error()  # Then convert to probability of error
        else:
            raise ValueError("No Phred score information set.")

    def _phred_char_to_q_score_int(self) -> int:
        if self.char is None:
            raise ValueError("Phred char is not set.")
        if self.q_score is None:
            self.q_score = ord(self.char.encode("utf-8")) - 33
        return self.q_score

    def _phred_q_score_int_to_prob_error(self) -> float:
        if self.q_score is None:
            raise ValueError("Phred Q-score is not set.")
        if self.prob_error is None:
            self.prob_error = 10 ** (-self.q_score / 10)
        return self.prob_error

    def _phred_q_score_int_to_char(self) -> str:
        if self.q_score is None:
            raise ValueError("Phred Q-score is not set.")
        if self.char is None:
            self.char = chr(self.q_score + 33).encode("ascii").decode()
        return self.char

    def _prob_error_to_phred_q_score_int(self) -> int:
        if self.prob_error is None:
            raise ValueError("Probability of error is not set.")
        if self.q_score is None:
            self.q_score = round(-10 * log10(self.prob_error))
        return self.q_score


class PhredString:
    """
    Class to handle strings of Phred quality scores.
    """
    def __init__(self, phred_string: str = None, phred_chars: List[NucleotideQuality] = None):
        self.phred_string = phred_string
        self.phred_chars = phred_chars
        if not any((
                (phred_string is not None),
                (phred_chars is not None),
        )):
            raise ValueError("At least one of the following must be set: phred_string, phred_chars.")
        if phred_chars is None:
            self.phred_chars = [NucleotideQuality(char=char) for char in phred_string]
        if phred_string is None:
            self.phred_string = "".join([pos_phred.to_phred_char() for pos_phred in self.phred_chars])

    def __len__(self):
        return len(self.phred_chars)

    def __getitem__(self, item):
        return self.phred_chars[item]

    def __str__(self):
        return self.phred_string

    def to_phred_string(self) -> str:
        return "".join([pos_phred.to_phred_char() for pos_phred in self.phred_chars])

    def to_phred_chars(self) -> List[str]:
        return [pos_phred.to_phred_char() for pos_phred in self.phred_chars]

    def to_q_scores(self) -> List[int]:
        return [pos_phred.to_q_score() for pos_phred in self.phred_chars]

    def to_prob_errors(self) -> List[float]:
        return [pos_phred.to_prob_error() for pos_phred in self.phred_chars]
    
    def adjust_values(self, adjustment_list: List[float]):
        assert len(adjustment_list) == len(self.phred_chars), "Adjustment list must be the same length as Phred string."
        for i, adjustment in enumerate(adjustment_list):
            self.phred_chars[i] = NucleotideQuality(
                prob_error=clamp((self.phred_chars[i].to_prob_error() * adjustment), 0, 1),
            )


def avg_phred_strings(phred_strings: List[PhredString]) -> PhredString:
    # We need to assert that all the strings are the same length
    assert len(set([len(phred_string) for phred_string in phred_strings])) == 1, \
        "Phred strings are not all the same length."
    phred_prob_errors = [phred.to_prob_errors() for phred in phred_strings]
    # Now we want to calculate the average probability of error for each position
    avg_prob_errors = []
    for nucl_pos in range(len(phred_prob_errors[0])):
        pos_prob_errors = [ind_phred_prob_errors[nucl_pos] for ind_phred_prob_errors in phred_prob_errors]
        avg_prob_errors.append(sum(pos_prob_errors) / len(pos_prob_errors))
    # Now we need to convert these back to Phred scores
    avg_phred_scores = [NucleotideQuality(prob_error=prob_error).to_phred_char()
                        for prob_error in avg_prob_errors]
    return PhredString("".join(avg_phred_scores))


def avg_raw_phred_strings(phred_strings: List[str]) -> PhredString:
    return avg_phred_strings([PhredString(phred) for phred in phred_strings])


def avg_raw_phred_strings_to_raw_string(phred_strings: List[str]) -> str:
    return avg_raw_phred_strings(phred_strings).to_phred_string()


def avg_phred_char_strings(phred_chars: List[str]) -> str:
    return avg_phred_strings([PhredString(phred_char) for phred_char in phred_chars]).to_phred_string()


def show_phred_table():
    """
    Show the Phred score table.
    """
    for val in range(0, 42+1):
        i = NucleotideQuality(q_score=val)
        print(f"Q:{i.to_q_score():>3}, ASCII: {val+33},"
              f"  Phred: {i.to_phred_char()}, Prob Error: {i.to_prob_error():>7.2g}")


def test_main(test_phreds, print_stuff=True):
    if print_stuff:
        show_phred_table()
    avg_phred = avg_raw_phred_strings(test_phreds)
    
    phred_strings = [PhredString(phred) for phred in test_phreds]
    
    spacer = 7
    if print_stuff:
        for phred_string in set(phred_strings):
            print(*[f"{i:^{spacer}}" for i in phred_string.to_phred_chars()])
            print(*[f"{i:^{spacer}.2g}" for i in phred_string.to_prob_errors()])
        print("-" * ((spacer + 1) * len(avg_phred)))
        print(*[f"{i:^{spacer}}" for i in avg_phred.to_phred_chars()])
        print(*[f"{i:^{spacer}.2g}" for i in avg_phred.to_prob_errors()])


if __name__ == '__main__':
    tester_phreds = ("8743233434/**'''',...442355*('&%%&&&('()./()))*+/))),--.874*)449:>>7('&''((*++-",
                     "8743233434/**'''',...442355*('&%%&&&('()./()))*+/))),--.874*)449:>>7('&''((*++-",
                     "/007(((((./,*&&'*+-('),++,3675441,+(%%))''(/.)(&'()*-334+**+*0//00/-.//3/-*&%%%",
                     "/007(((((./,*&&'*+-('),++,3675441,+(%%))''(/.)(&'()*-334+**+*0//00/-.//3/-*&%%%",
                     "/007(((((./,*&&'*+-('),++,3675441,+(%%))''(/.)(&'()*-334+**+*0//00/-.//3/-*&%%%")
    
    test_main(tester_phreds)
    
    # For very rough speed testing:
    for i in tqdm(range(100), desc=f"Testing Speed while averaging {len(tester_phreds):>4} Phred strings"):
        test_main(tester_phreds, print_stuff=False)
        
    # tester_phreds *= 10  # So we are now testing 50 phred strings per averaging step
    # for i in tqdm(range(100), desc=f"Testing Speed while averaging {len(tester_phreds):>4} Phred strings"):
    #     test_main(tester_phreds, print_stuff=False)
    #     
    # tester_phreds *= 10  # So we are now testing 500 phred strings per averaging step
    # for i in tqdm(range(100), desc=f"Testing Speed while averaging {len(tester_phreds):>4} Phred strings"):
    #     test_main(tester_phreds, print_stuff=False)
        
    # tester_phreds *= 10  # So we are now testing 5000 phred strings per averaging step
    # for i in tqdm(range(100), desc=f"Testing Speed while averaging {len(tester_phreds):>4} Phred strings"):
    #     test_main(tester_phreds, print_stuff=False)
    
    # Speed seems pretty linear!!
    
    # multiplying:
    test_confidence_scores = [clamp(random.random(), 0.5, 1.0) for _ in range(len(tester_phreds[0]))]
    # TODO: Ensure that we double check that we have reciprocal confidence scores when we do the real thing!
    #       This is because we want to lower the Phred score for positions with lower confidence scores.
    #       And since we are adjusting the probability of error, we need to bring those closer to 1!
    recip_confidence_scores = [1 / score for score in test_confidence_scores]
    test_phred_1 = PhredString(tester_phreds[0])
    
    spacer = 7
    print("\n\nTesting multiplication:")
    print("Original:")
    print(*[f"{i:^{spacer}}" for i in test_phred_1.to_phred_chars()])
    print(*[f"{i:^{spacer}.2g}" for i in test_phred_1.to_prob_errors()])
    print("-" * ((spacer + 1) * len(test_phred_1)))
    print("Confidence scores:")
    print(*[f"{i:^{spacer}.2g}" for i in test_confidence_scores])
    print("Reciprocal Confidence scores:")
    print(*[f"{i:^{spacer}.2g}" for i in recip_confidence_scores])
    print("-" * ((spacer + 1) * len(test_phred_1)))
    
    og_string = test_phred_1.to_phred_string()
    test_phred_1.adjust_values(recip_confidence_scores)
    new_string = test_phred_1.to_phred_string()
    print("Adjusted:")
    print(*[f"{i:^{spacer}}" for i in test_phred_1.to_phred_chars()])
    print(*[f"{i:^{spacer}.2g}" for i in test_phred_1.to_prob_errors()])
    print("-" * ((spacer + 1) * len(test_phred_1)))
    print(f"Original string: {og_string}")
    print(f"Adjusted string: {new_string}")
