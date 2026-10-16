import numpy as np

from gena.raw import packet, structure


def test_packets_2_pdf_decodes_big_endian_event_counts():
    packets = np.zeros(1, dtype=structure.packet_struct)
    packets["ep1"]["Evt_Num_Tot"] = 1
    packets["ep1"]["Evt_Num_Eff"] = 0x1234
    packets["ep1"]["Evt_Num_Aux"] = 0xFFFF

    pdf = packet.packets_2_pdf(packets)

    assert pdf[["Evt_Num_Tot", "Evt_Num_Eff", "Evt_Num_Aux"]].iloc[0].tolist() == [1, 0x1234, 0xFFFF]
