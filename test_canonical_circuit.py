import unittest

import hardware


class CanonicalCircuitTests(unittest.TestCase):
    def test_n6_fused_circuit_matches_reported_targets(self):
        n = 6
        steps = 100
        circuit = hardware.build_tfim_circuit_2nd_order(
            n, J=1.0, h=0.5, dt=0.025, steps=steps, periodic=True
        )
        counts = circuit.count_ops()
        gate_count = circuit.size()
        depth = circuit.depth()
        rx_layers = counts.get("rx", 0) // n

        targets = {"gates": 1218, "depth": 303, "rx_layers": steps + 1}
        actual = {"gates": gate_count, "depth": depth, "rx_layers": rx_layers}
        print(f"N=6 fused circuit actual={actual}; ops={dict(counts)}; targets={targets}")

        self.assertLess(abs(gate_count - targets["gates"]) / targets["gates"], 0.03)
        self.assertLess(abs(depth - targets["depth"]) / targets["depth"], 0.03)
        self.assertEqual(rx_layers, targets["rx_layers"])
        self.assertLess(gate_count, 1812)
        self.assertLess(depth, 402)


if __name__ == "__main__":
    unittest.main()
