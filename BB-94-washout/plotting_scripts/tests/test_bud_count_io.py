"""Checks that decoding cannot silently join the wrong samples or counts."""

from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bud_count_io import InputError, load_pairs
from decode_bud_counts import main as decode_main


class DecodingTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)

    def tearDown(self):
        self.temp.cleanup()

    def pair(self, records, counts, name="sample", delimiter="\t"):
        record = self.root / (name + "-record.txt")
        count = self.root / (name + "-count.txt")
        record.write_text("original_file_name" + delimiter + "converted_file_name\n" + "\n".join(delimiter.join(row) for row in records) + "\n")
        count.write_text("file_name" + delimiter + "bud_number\n" + "\n".join(delimiter.join(row) for row in counts) + "\n")
        return record, count

    def test_join_is_by_hash_and_subject_persists_across_dates(self):
        first = self.pair(
            [("210325-SMG-ctrl-t00_1.nd2", "hash1.nd2"), ("210325-SMG-ctrl-t00_2.nd2", "hash2.nd2")],
            [(" hash2.nd2 ", "7"), ("/images/hash1.nd2", "3")]
        )
        second = self.pair([("210326-SMG-ctrl-t22_1.nd2", "hash3.nd2")], [("hash3.nd2", "12")], name="later", delimiter=",")
        data = load_pairs([second, first], experiment="cohort")
        self.assertEqual(data.bud_number.tolist(), [3, 12, 7])
        self.assertEqual(data.time_h.tolist(), [0, 22, 0])
        self.assertEqual(data.subject_id.nunique(), 2)
        self.assertEqual(data.subject_id.iloc[0], data.subject_id.iloc[1])

    def test_missing_and_extra_counts_fail(self):
        pair = self.pair([("210325-SMG-ctrl-t00_1.nd2", "hash1.nd2")], [("wrong.nd2", "2")])
        with self.assertRaisesRegex(InputError, "missing counts.*counts without a record"):
            load_pairs([pair])

    def test_basename_collision_fails(self):
        pair = self.pair(
            [("210325-SMG-ctrl-t00_1.nd2", "hash1.nd2"), ("210325-SMG-ctrl-t00_2.nd2", "hash2.nd2")],
            [("a/hash1.nd2", "2"), ("b/hash1.nd2", "3")]
        )
        with self.assertRaisesRegex(InputError, "duplicate file_name"):
            load_pairs([pair])

    def test_invalid_counts_fail(self):
        for value in ("", "NaN", "inf", "-1", "2.5", "two", "9223372036854775808"):
            with self.subTest(value=value):
                pair = self.pair([("210325-SMG-ctrl-t00_1.nd2", "hash1.nd2")], [("hash1.nd2", value)])
                with self.assertRaises(InputError):
                    load_pairs([pair])

    def test_duplicate_gland_time_cells_fail_even_for_different_dates(self):
        pair = self.pair(
            [("210325-SMG-ctrl-t00_1.nd2", "hash1.nd2"), ("210326-SMG-ctrl-t00_1.nd2", "hash2.nd2")],
            [("hash1.nd2", "2"), ("hash2.nd2", "3")]
        )
        with self.assertRaisesRegex(InputError, "Duplicate experiment/group/gland/time"):
            load_pairs([pair])

    def test_metadata_overrides_nonstandard_names_and_is_exact(self):
        pair = self.pair([("sample.tif", "random.tif")], [("random.tif", "2")])
        metadata = self.root / "metadata.csv"
        metadata.write_text("original_file_name,experiment,group,gland,time_h\nsample.tif,cohort2,BB94,03,22\n")
        row = load_pairs([pair], metadata_path=metadata).iloc[0]
        self.assertEqual((row.experiment, row.group, row.gland, row.time_h), ("cohort2", "BB94", "03", 22))
        metadata.write_text(metadata.read_text() + "unused.tif,cohort2,BB94,04,22\n")
        with self.assertRaisesRegex(InputError, "unused mappings"):
            load_pairs([pair], metadata_path=metadata)

    def test_repeated_pair_fails(self):
        pair = self.pair([("210325-SMG-ctrl-t00_1.nd2", "hash.nd2")], [("hash.nd2", "2")])
        with self.assertRaisesRegex(InputError, "Repeated original_file_name"):
            load_pairs([pair, pair])

    def test_days_convert_to_hours_and_control_alias_is_canonical(self):
        pair = self.pair(
            [("210401-plate2-control-day0_5.nd2", "hash1.nd2"),
             ("20210402-PLATE2-CONTROL-DAY1.5_5.nd2", "hash2.nd2")],
            [("hash2.nd2", "12"), ("hash1.nd2", "3")]
        )
        data = load_pairs([pair])
        self.assertEqual(data.group.tolist(), ["ctrl", "ctrl"])
        self.assertEqual(data.gland.tolist(), ["plate2_5", "plate2_5"])
        self.assertEqual(data.plate.tolist(), ["plate2", "plate2"])
        self.assertEqual(data.time_unit.tolist(), ["day", "day"])
        self.assertEqual(data.time_value.tolist(), [0, 1.5])
        self.assertEqual(data.time_h.tolist(), [0, 36])
        self.assertEqual(data.subject_id.nunique(), 1)

    def test_same_local_gland_on_different_plates_is_distinct(self):
        pair = self.pair(
            [("210401-plate1-control-day0_5.nd2", "hash1.nd2"),
             ("210401-plate2-control-day0_5.nd2", "hash2.nd2")],
            [("hash1.nd2", "3"), ("hash2.nd2", "4")]
        )
        data = load_pairs([pair])
        self.assertEqual(data.gland.tolist(), ["plate1_5", "plate2_5"])
        self.assertEqual(data.subject_id.nunique(), 2)

    def test_hour_convention_keeps_gland_and_subject_identity(self):
        pair = self.pair(
            [("210325-SMG-BB94-t22_8.nd2", "hash1.nd2")], [("hash1.nd2", "5")]
        )
        row = load_pairs([pair], experiment="washout").iloc[0]
        self.assertEqual(row.gland, "8")
        self.assertEqual(row.subject_id, "washout|BB94|8")
        self.assertEqual(row.plate, "")
        self.assertEqual((row.time_unit, row.time_value, row.time_h), ("hour", 22, 22))

    def test_mixed_structural_conventions_decimal_units_and_bb94_alias(self):
        pair = self.pair(
            [("20210401-smg-bb-94-Day0.5_3.tif", "hash1.tif"),
             ("20210401-PlAtE3-bb94-T2.5_3.tif", "hash2.tif")],
            [("hash1.tif", "3"), ("hash2.tif", "4")]
        )
        data = load_pairs([pair]).set_index("gland")
        self.assertEqual(data.loc["3", "group"], "BB94")
        self.assertEqual(data.loc["3", "time_h"], 12)
        self.assertEqual(data.loc["3", "time_unit"], "day")
        self.assertEqual(data.loc["plate3_3", "time_h"], 2.5)
        self.assertEqual(data.loc["plate3_3", "time_unit"], "hour")

    def test_negative_or_malformed_day_is_rejected(self):
        for time in ("-1", "NaN", "1.2.3"):
            with self.subTest(time=time):
                pair = self.pair([("210401-plate2-control-day{}_5.nd2".format(time), "hash.nd2")], [("hash.nd2", "3")])
                with self.assertRaisesRegex(InputError, "Cannot parse"):
                    load_pairs([pair])

    def test_aliases_cannot_create_duplicate_gland_time_observations(self):
        pair = self.pair(
            [("210401-plate2-control-day1_5.nd2", "hash1.nd2"),
             ("210402-PLATE02-ctrl-t24_5.nd2", "hash2.nd2")],
            [("hash1.nd2", "3"), ("hash2.nd2", "4")]
        )
        with self.assertRaisesRegex(InputError, "Duplicate experiment/group/gland/time"):
            load_pairs([pair])

    def test_metadata_time_and_labels_override_filename_units(self):
        pair = self.pair([("210401-plate2-control-day1_5.nd2", "hash.nd2")], [("hash.nd2", "3")])
        metadata = self.root / "metadata.csv"
        metadata.write_text("original_file_name,group,gland,time_h,time_unit,time_value,plate\n"
                            "210401-plate2-control-day1_5.nd2,Control custom,plate2-local5,22,day,1,ignored\n")
        row = load_pairs([pair], metadata_path=metadata).iloc[0]
        self.assertEqual((row.group, row.gland), ("Control custom", "plate2-local5"))
        self.assertEqual((row.time_h, row.time_unit, row.time_value, row.plate), (22, "hour", 22, ""))

    def test_d_and_day_spellings_have_same_units_and_gland_identity(self):
        pair = self.pair(
            [("210401-plate2-control-d0_5.nd2", "hash1.nd2"),
             ("210402-plate2-control-day1_5.nd2", "hash2.nd2"),
             ("210403-plate2-control-D2_5.nd2", "hash3.nd2")],
            [("hash1.nd2", "3"), ("hash2.nd2", "4"), ("hash3.nd2", "12")]
        )
        data = load_pairs([pair])
        self.assertEqual(data.time_h.tolist(), [0, 24, 48])
        self.assertEqual(data.time_value.tolist(), [0, 1, 2])
        self.assertEqual(data.time_unit.tolist(), ["day"] * 3)
        self.assertEqual(data.subject_id.nunique(), 1)

    def test_day_alias_cannot_duplicate_same_gland_observation(self):
        pair = self.pair(
            [("210401-plate2-control-day1_5.nd2", "hash1.nd2"),
             ("210401-plate2-control-d1_5.nd2", "hash2.nd2")],
            [("hash1.nd2", "3"), ("hash2.nd2", "4")]
        )
        with self.assertRaisesRegex(InputError, "Duplicate experiment/group/gland/time"):
            load_pairs([pair])

    def test_decoder_refuses_existing_output_and_input_overwrite(self):
        pair = self.pair([("210325-SMG-ctrl-t00_1.nd2", "hash.nd2")], [("hash.nd2", "2")])
        output = self.root / "decoded.csv"
        args = ["--pair", str(pair[0]), str(pair[1]), "--output", str(output)]
        self.assertEqual(decode_main(args), 0)
        with self.assertRaises(SystemExit):
            decode_main(args)
        with self.assertRaises(SystemExit):
            decode_main(["--pair", str(pair[0]), str(pair[1]), "--output", str(pair[0]), "--overwrite"])


if __name__ == "__main__":
    unittest.main()
