from functions import *
import streamlit as st


def add_to_history(G, base, y, decoded_y):
    st.session_state.history.append(
        {"G": G, "base": base, "y": y, "decoded_y": decoded_y})


def rerender_history():
    with st.expander("History"):
        for entry in st.session_state.history:
            text = ""
            for d in entry:
                if d == "G":
                    text += "G:\n"
                    for l in entry[d]:
                        text += ""
                        for t in l:
                            text += f"{str(t)} "
                        text += "\n"
                else:
                    text += f"{d}:\t{entry[d]}\n"
            st.code(text)


if __name__ == "__main__":
    st.write("""
        # Linear Code Tools
        """)

    G = """1 0 0 0 1 1
0 1 0 1 0 1
0 0 1 1 1 0"""

    if 'history' not in st.session_state:
        st.session_state.history = []

    with st.form("Generator form"):
        G = st.text_area("Input your Generator matrix(G):",
                         value=G, height=150)
        base = st.text_input("Base:", value="2")
        y = st.text_input("Input your recieved data(y):", value="1 1 1 1 0 1")
        submit_button = st.form_submit_button(label='Submit')

    if submit_button:
        G = [[int(g) for g in l.replace(" ", "")] for l in G.split("\n")]
        base = int(base)
        y = [int(x) for x in y.replace(" ", "")]

        standard_form_G = transform_to_standard_form(G, base)
        G2 = np.array(standard_form_G)
        k, n = G2.shape

        if len(y) != n:
            st.warning("wrong recieved code")
            decoded_y = "none"

        C = generate_codewords(G2, base)
        st.code(f"""Codewords (C):\n{C}\n\nminimum distance:\t
                {calculate_minimum_distance(C)}""")

        H = parity_check_matrix_generator(G2, base)
        text = "Parity Check matrix (H):\n"
        for l in H:
            text += ""
            for t in l:
                text += f"{str(t)} "
            text += "\n"
        st.code(text)

        syndrome_lookup_table = {}
        syndromes_cosets = generate_syndrome_table(H, base)
        text = "syndrome\t\tcosets\n"
        for syndrome in syndromes_cosets:
            text += f"{syndrome}\t\t"
            all_cosets_of_syndrome = syndromes_cosets[syndrome]
            text += f"{[list(x) for x in all_cosets_of_syndrome]}\n"
            syndrome_lookup_table[syndrome] = arrays_with_lowest_sum(
                all_cosets_of_syndrome)
        st.code(text)

        text = "Syndrome lookup table:\nsyndrome z\t\t\tcoset leader f(z)\n"
        for syndrome in syndrome_lookup_table:
            text += f"{syndrome}\t\t\t{syndrome_lookup_table[syndrome]}\n"
        st.code(text)

        try:
            z = calculate_syndrome(y, H, base)
            fz = syndrome_lookup_table[tuple(z)]
            if len(fz) >= 2:
                st.warning("there are more than one error")
                add_to_history(G, base, y, "none")
            else:
                decoded_y = [(((p-q)+base) % base) for p, q in zip(y, fz[0])]
                if y == decoded_y:
                    st.code(f"Code was correct:\n{decoded_y}")
                else:
                    st.code(f"Correct code:\n{decoded_y}")
                add_to_history(G, base, y, decoded_y)
        except:
            add_to_history(G, base, y, "none")
        rerender_history()
