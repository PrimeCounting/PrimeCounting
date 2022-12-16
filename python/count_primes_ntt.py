import time
import array
import itertools
import math
import sympy


def get_primes(max_n):
    """Find primes < `max_n` using a sieve."""
    sieve = array.array('b', [1]) * max_n
    sieve[0] = sieve[1] = 0
    result = []
    for i in range(max_n):
        if sieve[i] == 1:
            # `i` is a prime.
            result.append(i)
            for j in range(i * i, max_n, i):
                sieve[j] = 0
    return result


class CountPrimes:
    """Count primes < `max_n`."""

    def __init__(self, max_n, profile=True):
        self.max_n = max_n

        # Track the current time for profiling.
        self.profile = profile
        self.previous_time = time.time()

        self.max_exact_prime = int(self.max_n ** 0.5) + 1
        # Assert that there are no rounding issues.
        # This assertion is critical.
        assert self.max_exact_prime ** 2 > self.max_n
        # This assertion not important, we just need self.max_exact_prime <= self.max_n.
        assert (self.max_exact_prime - 1) ** 2 <= self.max_n

        # Primes that are required for mobius inclusion-exclusion.
        self.exact_primes = get_primes(self.max_exact_prime)

        # The number `n` will be projected to `floor(ln(n) * self.log_multiplier)`.
        # So the logarithmic precision is `1 / self.log_multiplier`.
        # `log_multiplier` is chosen to balance the approximation and correction times.
        # The constant here was chosen heuristically.
        log_max_n = math.log(self.max_n)
        self.log_multiplier = 0.7 * self.max_exact_prime / log_max_n

        # Number of interesting log bins.
        # Bins beyond this will be discarded.
        self.num_log_bins = self.log_bin(max_n) + 1

        # Maximal exponent of a prime that fits in the log bins.
        # Note that in practice non-squarefree numbers are ignored by mobius so instead of `max_prime_exponent` we are
        # interested in the maximal number of unique prime factors a number may have.
        # TODO: Replace this bound by the stricter value to avoid unneeded padding.
        self.max_prime_exponent = (self.num_log_bins - 1) // self.log_bin(2) + 1

        # Find the maximal value that appears during the computation of the mobius function.
        # TODO: Partition the primes according to size to reduce `max_value` significantly
        max_value = self.max_exact_prime ** self.max_prime_exponent
        self.num_ntt_bins = 2 ** (self.log_bin(max_value).bit_length() + 1)

        # The modulus of our NTT/FFT.
        self.ntt_prime = (3 << 30) + 1
        assert pow(7, self.ntt_prime - 1, self.ntt_prime) == 1, "Should be prime."

        # Assert that p-1 is divisible by enough 2 factors.
        ntt_prime_totient = self.ntt_prime - 1
        ntt_prime_totient_two_factors = ntt_prime_totient & ~(ntt_prime_totient - 1)
        assert self.num_ntt_bins <= ntt_prime_totient_two_factors

        # Assert that `ntt_prime` is large enough to contain the result.
        # Note that we could use the Riemann hypothesis to bound the result in a small interval so the modulus could be
        # much smaller: sqrt(max_n) * log(max_n).
        approximate_num_primes = 1.01 * self.max_n / math.log(self.max_n) + 100
        assert approximate_num_primes < self.ntt_prime

        self.timing("Constructor")

    def timing(self, title):
        if not self.profile:
            return
        current_time = time.time()
        print(f"{title} took {current_time - self.previous_time:.4} seconds.")
        self.previous_time = current_time

    def log_bin(self, n):
        return int(math.log(n) * self.log_multiplier)

    def max_log_error(self):
        """Error is one-directional towards the zero."""
        return 1 / self.log_multiplier

    def new_counter(self):
        """An empty array of log bins."""
        return [0] * self.num_ntt_bins

    def convolve_update(self, side0, side1, coeff, result):
        """Convolve two logarithmic counters and accumulate the result."""
        for n0, c0 in enumerate(side0):
            for n1, c1 in enumerate(side1):
                n_new = n0 + n1
                if n_new < self.num_log_bins:
                    result[n_new] += coeff * c0 * c1

    def convolve(self, side0, side1, coeff):
        result = self.new_counter()
        self.convolve_update(side0, side1, coeff, result)
        return result

    def count_primes_approximate(self):
        primes_counter = self.new_counter()
        for p in self.exact_primes:
            primes_counter[self.log_bin(p)] += 1
        self.timing("Exact prime counting")

        primes_ntt = sympy.discrete.transforms.ntt(primes_counter, prime=self.ntt_prime)
        self.timing("Exact prime NTT")

        # The NTT of the counter of mobius(n) for `max_exact_prime`-smooth numbers.
        # We call this function the smooth mobius.
        # We compute it by accumulating products of exact primes with alternating signs.
        mobius_ntt = []

        inverse = [pow(x, self.ntt_prime - 2, self.ntt_prime) for x in range(self.max_prime_exponent)]
        for ntt_bin in range(self.num_ntt_bins):
            # prime_product_ntts[multiplicity] is the NTT coefficient of products of `multiplicity` unique primes.
            # Initialize with `1`, the empty product. All its NTT coefficients equals to 1.
            prime_product_ntts = [1]

            # We use the recursion (aka Newton's Identities):
            # n * prime_product_counters[n] =
            #   conv(prime_power_counters[1], prime_product_counters[n-1]) -
            #   conv(prime_power_counters[2], prime_product_counters[n-2]) +
            #   conv(prime_power_counters[3], prime_product_counters[n-3]) -
            #   ...
            for multiplicity in range(1, self.max_prime_exponent):
                s = 0
                # Inclusion-Exclusion a la Newton.
                # Iterate in reverse because exponent = 1 should have a plus sign.
                for exponent in reversed(range(1, multiplicity + 1)):
                    # NTT of the counter of prime**exponent is simply a shuffled version of `primes_ntt`.
                    # The `ntt_bin` is multiplied by the exponent, modulo the number of NTT bins.
                    prime_exponent_ntt_bin = exponent * ntt_bin % self.num_ntt_bins
                    s = primes_ntt[prime_exponent_ntt_bin] * prime_product_ntts[-exponent] - s
                # Normalize by the number of symmetries.
                prime_product_ntts.append((s * inverse[multiplicity]) % self.ntt_prime)

            # Compute mobius by accumulating products of exact primes with alternating signs.
            mobius_accumulator = 0
            # Iterate in reverse because prime_product_ntts[0] should have a plus sign.
            for prime_product_ntt in reversed(prime_product_ntts):
                mobius_accumulator = prime_product_ntt - mobius_accumulator
            mobius_ntt.append(mobius_accumulator)
        self.timing("Newton's Identities + Mobius")

        # A counter of all numbers.
        all_numbers = self.new_counter()
        n_previous = 0
        for log_bin in range(self.num_log_bins):
            # Find the maximal n that is contained in the current log bin.
            # Start by jumping the logarithmic interval of a single bin.
            n = int(n_previous * math.exp(1 / self.log_multiplier))
            if n == 0:
                n += 1
            assert self.log_bin(n) <= log_bin
            while self.log_bin(n + 1) == log_bin:
                n += 1
            all_numbers[log_bin] = n - n_previous
            n_previous = n
        self.timing("All numbers counting")

        all_numbers_ntt = sympy.discrete.transforms.ntt(all_numbers, prime=self.ntt_prime)
        self.timing("All numbers NTT")

        # Convolving mobius and all numbers removes numbers that are divisible by any of the exact primes.
        # Since we took exact primes up to sqrt, this leaves only prime numbers larger than the exact primes.
        large_primes_ntt = [a * b for a, b in zip(mobius_ntt, all_numbers_ntt)]
        self.timing("Final convolution")

        large_primes = sympy.discrete.transforms.intt(large_primes_ntt, prime=self.ntt_prime)
        num_primes = sum(large_primes[:self.num_log_bins])  # Counts primes above self.max_exact_prime but also 1.
        num_primes %= self.ntt_prime
        num_primes -= 1
        num_primes += len(self.exact_primes)
        self.timing("Inverse NTT")
        return num_primes

    def count_corrections_for_number(self, n, exact_prime_divisors):
        """
        Cancel any contribution of the number `n`.
        `exact_prime_divisors` contains `exact_primes` that divides `n`, no multiplicity.
        """
        # We assume `n` to be outside of the target range, so any contribution it has made should be undone.
        assert n >= self.max_n
        if not exact_prime_divisors:
            # `n` does not have any exact prime divisor.
            # `n` was only counted directly in `all_numbers` if it's in the log bins.
            return self.log_bin(n) < self.num_log_bins

        # Contribution made to log bin < log_bin_threshold should be counted.
        # We will modify this threshold along the computation.
        log_bin_threshold = self.num_log_bins

        exact_prime_product = 1
        for p in exact_prime_divisors:
            log_bin_threshold -= self.log_bin(p)
            exact_prime_product *= p
        if log_bin_threshold <= 0:
            return 0

        # Accumulates the contributions of the number `n`.
        correction = 0

        # Find all squarefree and smooth divisors of `n` by iteratively adding each prime.
        # We only look for squarefree and smooth because these are the only numbers included by our smooth-mobius
        # function.
        # The following three arrays are aligned so their i-th element describes the same divisor.
        n_divisors = []
        n_divisors_log_bin_threshold = []
        n_divisors_mobius = []

        def push_divisor(divisor, divisor_log_bin_threshold, divisor_mobius):
            # Check whether `n = divisor * leftover` was counted correctly where divisor
            # was in the `mobius` counter and `leftover` was in the `all_number` counter.
            leftover = n // divisor
            log_bin_leftover = self.log_bin(leftover)
            if log_bin_leftover >= divisor_log_bin_threshold:
                return 0
            n_divisors.append(divisor)
            n_divisors_log_bin_threshold.append(divisor_log_bin_threshold)
            n_divisors_mobius.append(divisor_mobius)
            # The number `n` when expressed as `n = divisor * leftover` counted with sign `divisor_mobius`.
            # We assumed the leftover was small enough to be included in `all_numbers`.
            # That may not be the case if `max_n` is very small.
            assert log_bin_leftover < self.num_log_bins
            return divisor_mobius

        # Initialization.
        correction += push_divisor(exact_prime_product, log_bin_threshold, (-1) ** len(exact_prime_divisors))

        # This costs 2**len(exact_prime_divisors).
        # TODO: By shifting the value of each prime by its rounded log, the task is reduced to counting the number of
        # subsets of [log(p) - rounded_log(p) for p in exact_prime_divisors] whose sum is smaller than a threshold.
        # This direction enables two improvements:
        # 1. Backtrack when computing the subsets so only good subsets (with a bounded sum) are generated.
        # 2. Split the primes into two sets and use MITM to count good subsets.
        # Still need to take care of signs (mobius) and verify that it works.
        for p in exact_prime_divisors:
            assert n % p == 0
            log_bin_p = self.log_bin(p)
            for i in range(len(n_divisors)):
                correction += push_divisor(
                    divisor=n_divisors[i] // p,
                    divisor_log_bin_threshold=n_divisors_log_bin_threshold[i] + log_bin_p,
                    divisor_mobius=-n_divisors_mobius[i],
                )
        return correction

    def count_corrections(self):
        max_cumulative_log_error = self.max_log_error() * self.max_prime_exponent
        # Numbers close to `max_n` up `max_cumulative_log_error` may have been counted incorrectly.
        min_n_suspected = self.max_n
        max_n_suspected = int(self.max_n * math.exp(max_cumulative_log_error)) + 2

        # Sieve to find exact prime divisors in the suspected range.
        # Note that we only sieve `exact_primes`. Suspected numbers are larger than `max_n` so they may be a product of
        # two primes larger than the exact primes.
        # This is not a problem since we do not wish to factorize the numbers, we just need divisors that were included
        # in the smooth mobius, meaning that they are `exact_prime`-smooth.
        prime_divisors_sieve = [[] for _ in range(min_n_suspected, max_n_suspected)]
        for p in self.exact_primes:
            start = ((-min_n_suspected) % p)
            for offset in range(start, max_n_suspected - min_n_suspected, p):
                prime_divisors_sieve[offset].append(p)
        self.timing("Sieving")

        total_correction = sum(
            self.count_corrections_for_number(n, prime_divisors_sieve[n - min_n_suspected])
            for n in range(min_n_suspected, max_n_suspected)
        )
        self.timing("Corrections")
        return total_correction

    def count_primes(self):
        return self.count_primes_approximate() - self.count_corrections()


def compare(max_n):
    tt = time.time()
    result = CountPrimes(max_n).count_primes()
    fast_time = time.time() - tt
    tt = time.time()
    reference = len(get_primes(max_n))
    naive_time = time.time() - tt
    print(
        f"Count primes up to {max_n}: fast={result} vs naive={reference}. "
        f"Took {fast_time:.3} vs {naive_time:.3} seconds.")
    assert result == reference
    return result


def run(max_n):
    tt = time.time()
    result = CountPrimes(max_n).count_primes()
    fast_time = time.time() - tt
    print(f"Count primes up to {max_n}: {result = }. Took {fast_time:.3} seconds.")
    return result


if __name__ == "__main__":
    # for max_n in range(2, 10**3):
    # for max_n in range(40, 10 ** 5, 177):
    #     compare(max_n)
    compare(10 ** 8)
