import numpy as np
from pymatgen.analysis.magnetism import CollinearMagneticStructureAnalyzer

from atomate.utils.utils import get_logger

logger = get_logger(__name__)


def procure_response_dict(
    struct_final,
    num_perturb_sites,
    incar_dict,
    outcar_dict,
    inv_block_dict,
    response_dict,
    perturb_dict,
    rkey,
    # keys,
    ldaul_vals,
    analyzer_gs,
    calcs_skipped,
):
    """
    Function to gather response data, in preparation for linear regression.
    This data is organized into `response_dict`.
    """

    # perform magnetic ordering analysis
    analyzer_output = CollinearMagneticStructureAnalyzer(struct_final, threshold=0.61)
    magnet_order = analyzer_output.ordering.value
    # if rkey == keys[0]:  # store ground state ordering
    #     magnet_order_gs = magnet_order

    # check if ordering matches ground state configuration
    if analyzer_gs:
        if not analyzer_gs.matches_ordering(struct_final):
            # use_calc = False
            calcs_skipped.append(
                {
                    "ICHARG": incar_dict.get("ICHARG", 0),
                    "ISPIN": incar_dict.get("ISPIN", 1),
                    "LDAUU": incar_dict["LDAUU"].copy(),
                    "LDAUJ": incar_dict["LDAUJ"].copy(),
                }
            )

    for i in range(num_perturb_sites):
        specie = struct_final[i].specie
        ldaul = ldaul_vals[i]

        orbital = inv_block_dict[str(ldaul)]
        perturb_dict.update({f"site{i}": {"specie": str(specie), "orbital": orbital}})

        # Obtain occupancy values
        n_tot = float(outcar_dict["charge"][i][orbital])
        # FIXME: Adapt for noncollinear
        m_z = float(outcar_dict["magnetization"][i][orbital])
        n_up = 0.5 * (n_tot + m_z)
        n_dn = 0.5 * (n_tot - m_z)

        v_up = float(incar_dict["LDAUU"][i])
        v_dn = float(incar_dict["LDAUJ"][i])

        response_dict[rkey][f"site{i}"]["Nup"].append(n_up)
        response_dict[rkey][f"site{i}"]["Ndn"].append(n_dn)
        response_dict[rkey][f"site{i}"]["Ntot"].append(n_tot)
        response_dict[rkey][f"site{i}"]["Mz"].append(m_z)

        response_dict[rkey][f"site{i}"]["Vup"].append(v_up)
        response_dict[rkey][f"site{i}"]["Vdn"].append(v_dn)

        response_dict[rkey]["magnetic order"].append(magnet_order)


def response_fit(x, y):
    """
    Function for fitting to response data. Returns: slope and associated error
    """

    (p, pcov) = np.polyfit(x, y, 1, cov=True)
    perr = np.sqrt(np.diag(pcov))

    return p, perr


def response_fit_stepped(x, y, tol=1.0e-6):
    """
    Function for fitting to response data
    - includes the "slope ~ zero" case for stepped data due to low precision
    Returns: slope and associated error
    """

    is_stepped = False
    step_id = -1

    y_sort = [y for y, _ in sorted(zip(y, x))]
    x_sort = [x for _, x in sorted(zip(y, x))]

    buff_size = 1  # must be gte three for first-order fit
    for i in range(buff_size, len(y_sort) - buff_size + 1):
        if np.std(y_sort[0:i]) < tol and np.std(y_sort[i:]) < tol:
            is_stepped = True
            step_id = i
            break

    if is_stepped:
        buff_max = 3  # must be >= three for first-order fit
        if step_id < buff_max:
            (p, perr) = response_fit(x[step_id:], y[step_id:])
        elif step_id > (len(y) - buff_max):
            (p, perr) = response_fit(x[0:step_id], y[0:step_id])
        else:
            (p1, p1err) = response_fit(x_sort[0:step_id], y_sort[0:step_id])
            (p2, p2err) = response_fit(x_sort[step_id:], y_sort[step_id:])
            p = 0.5 * (np.array(p1) + np.array(p2))
            perr = np.sqrt(0.5 * (np.array(p1err) ** 2 + np.array(p2err) ** 2))
        return p, perr
    else:
        (p, perr) = response_fit(x, y)
        return p, perr


def obtain_response_matrices(
    n_response,
    spin_polarized,
    response_dict,
    keys,
):
    """
    Function to compute self-consistent (SCF) and non-self-consistent (NSCF)
    linear response "chi" matrices; In addition to using linear regression
    to compute slopes about zero potential, the uncertainty associated
    with these values are also stored for subsequent error quantification.
    Returns: chi_matrix_nscf, chi_matrix_scf, chi_nscf_err, chi_scf_err
    """

    # Matrices for self-consistent and non-self-consistent responses
    # & associated element-wise errors
    chi_matrix_nscf = np.zeros([n_response, n_response])
    chi_matrix_scf = np.zeros([n_response, n_response])
    chi_nscf_err = np.zeros([n_response, n_response])
    chi_scf_err = np.zeros([n_response, n_response])

    # Compute response matrices using fitting function
    for ii in range(n_response):
        for jj in range(n_response):
            if spin_polarized:
                i, j = ii // 2, jj // 2
                si, sj = (
                    "up" if np.mod(ii, 2) == 0 else "dn",
                    "up" if np.mod(jj, 2) == 0 else "dn",
                )
                v_key = "V" + sj
                n_key = "N" + si
            else:
                i, j = ii, jj
                v_key = "Vup"
                n_key = "Ntot"

            if (
                response_dict[keys[1]][f"site{i}"][n_key]
                and response_dict[keys[2]][f"site{i}"][n_key]
                and response_dict[keys[1]][f"site{j}"][v_key]
                and response_dict[keys[2]][f"site{j}"][v_key]
            ):

                # gather NSCF & SCF response data
                v_nscf, n_nscf = [], []
                v_scf, n_scf = [], []
                for ll in [1, 2]:
                    for idx in range(len(response_dict[keys[ll]][f"site{j}"][v_key])):

                        v = response_dict[keys[ll]][f"site{j}"][v_key][idx]
                        n = response_dict[keys[ll]][f"site{i}"][n_key][idx]
                        # order = response_dict[keys[ll]]["magnetic order"][l]

                        # if order == magnet_order_gs:

                        isolated_response = v != 0.0

                        if isolated_response:
                            if ll == 1:
                                v_nscf.append(v)
                                n_nscf.append(n)
                            elif ll == 2:
                                v_scf.append(v)
                                n_scf.append(n)

                # Add ground state
                if (
                    response_dict[keys[0]][f"site{j}"][v_key]
                    and response_dict[keys[0]][f"site{i}"][n_key]
                ):
                    v = response_dict[keys[0]][f"site{j}"][v_key][0]
                    n = response_dict[keys[0]][f"site{i}"][n_key][0]
                    v_nscf.append(v)
                    n_nscf.append(n)
                    v_scf.append(v)
                    n_scf.append(n)

                try:
                    fit_nscf = response_fit(v_nscf, n_nscf)
                    chi_nscf, err_chi_nscf = fit_nscf[0][-2], fit_nscf[1][-2]
                    fit_scf = response_fit(v_scf, n_scf)
                    chi_scf, err_chi_scf = fit_scf[0][-2], fit_scf[1][-2]
                except Exception as exc:
                    chi_nscf, err_chi_nscf = float("nan"), float("nan")
                    chi_scf, err_chi_scf = float("nan"), float("nan")
                    logger.warning("Slope fitting fail", exc)

            else:
                chi_nscf, err_chi_nscf = float("nan"), float("nan")
                chi_scf, err_chi_scf = float("nan"), float("nan")

            chi_matrix_nscf[ii, jj] = chi_nscf
            chi_matrix_scf[ii, jj] = chi_scf
            chi_nscf_err[ii, jj] = err_chi_nscf
            chi_scf_err[ii, jj] = err_chi_scf

    return chi_matrix_nscf, chi_matrix_scf, chi_nscf_err, chi_scf_err


def inverse_matrix_uncertainty(matrix, matrix_covar):
    """
    Function to compute the element-wise error propagation in matrix inversion
    """

    m, n = matrix.shape
    if m != n:
        logger.warning("Matrix dimension error")
        return float("nan") * matrix, float("nan") * matrix

    matrixinv = np.linalg.inv(matrix)
    matrixinv_var = np.zeros([m, n])

    if m == 1 and n == 1:
        jacobian = -1 / matrix[0, 0] ** 2
        jacobians = [[jacobian]]
        sigma_f = jacobian * matrix_covar[0, 0] * jacobian
        matrixinv_var[0, 0] = sigma_f

        return matrixinv, matrixinv_var, jacobians

    # Function to determine the symbolic partial derivative of the
    # determinant w.r.t. matrix element
    def det_deriv(matrix, i, j):
        mij = np.delete(np.delete(matrix, i, 0), j, 1)
        partial = (-1) ** (i + j) * np.linalg.det(mij)
        return partial

    # Jacobians of each element of matrix inversion w.r.t.
    # original matrix elements
    jacobians = [[] for i in range(m)]

    det = np.linalg.det(matrix)
    for i in range(m):
        for j in range(n):
            mji = np.delete(np.delete(matrix, j, 0), i, 1)
            minor = (-1) ** (i + j) * np.linalg.det(mji)

            j_matrix = np.zeros([m, n])
            for k in range(m):
                for l in range(n):
                    det_p = det_deriv(matrix, k, l)

                    if k == j or l == i:
                        minor_p = 0.0
                    else:
                        kk, ll = k - 1 if k > j else k, l - 1 if l > i else l
                        minor_p = (-1) ** (i + j) * det_deriv(mji, kk, ll)

                    j_matrix[k, l] = (minor_p * det - minor * det_p) / det ** 2

            jacobians[i].append(j_matrix)

            j_vec = np.reshape(j_matrix, [m * n, 1])
            sigma_f = np.sum(np.dot(np.transpose(j_vec), np.dot(matrix_covar, j_vec)))
            matrixinv_var[i, j] = sigma_f

    return matrixinv, matrixinv_var, jacobians


def chi_inverse(chi, chi_err, method="full"):
    """
    Function to compute inverse of response matrix and associated
    element-wise uncertainty for point-wise, atom-wise,
    and full matrix inversion
    """

    n_response = len(chi)

    chi_block = chi.copy()
    chi_err_block = chi_err.copy()

    if method == "point":
        # diagonal 1x1
        for ii in range(n_response):
            for jj in range(n_response):
                if ii != jj:
                    chi_block[ii, jj], chi_err_block[ii, jj] = 0.0, 0.0
    elif method == "atom":
        # 2x2 block diagonal
        for ii in range(n_response):
            for jj in range(n_response):
                i, j = ii // 2, jj // 2
                if i != j:
                    chi_block[ii, jj], chi_err_block[ii, jj] = 0.0, 0.0
    elif method != "full":
        raise ValueError(
            "Unsupported method, method must be point (diagonal 1x1 inversion), "
            "atom (block 2x2 inverse), or full (full inverse)"
        )

    # Assume cross-covariances are zero
    chi_covar = np.diag(np.reshape(chi_err_block ** 2, [n_response * n_response]))

    (chi_inv, chi_inv_var, chi_inv_jacobs) = inverse_matrix_uncertainty(
        chi_block, chi_covar
    )

    return chi_block, chi_inv, chi_inv_var, chi_inv_jacobs


def compute_u_pointwise(
    site_index,
    f_matrix,
    f_matrix_err,
):
    """
    Function to compute Hubbard U value using point-wise (diagonal) inversion,
    in addition to the associated uncertainty value
    - based on the study by Linscott et. al.
    """

    i = site_index

    umat = f_matrix[2 * i : 2 * (i + 1), 2 * i : 2 * (i + 1)]
    umat_err = f_matrix_err[2 * i : 2 * (i + 1), 2 * i : 2 * (i + 1)]
    uval = 0.5 * np.sum(np.diag(umat))
    uval_err = 0.5 * np.sqrt(np.sum(np.diag(umat_err) ** 2))

    return uval, uval_err


def compute_uj_simple_two_by_two(
    site_index,
    f_matrix,
    f_matrix_err,
):
    """
    Function to compute Hubbard U and Hund J values using simple 2x2 formula,
    in addition to the associated uncertainty values
    - based on the study by Linscott et. al.
    """

    i = site_index

    umat = f_matrix[2 * site_index : 2 * (i + 1), 2 * i : 2 * (i + 1)]
    umat_err = f_matrix_err[2 * i : 2 * (i + 1), 2 * i : 2 * (i + 1)]

    uval = 0.25 * np.sum(umat)
    uval_err = 0.25 * np.sqrt(np.sum(umat_err ** 2))

    jmat = np.array([[-1, 1], [1, -1]]) * umat.copy()
    jmat_err = umat_err.copy()
    jval = 0.25 * np.sum(jmat)
    jval_err = 0.25 * np.sqrt(np.sum(jmat_err ** 2))

    return uval, uval_err, jval, jval_err


def compute_uj_scaled_two_by_two(
    site_index,
    f_matrix,
    f_matrix_err,
    chi_matrix_scf,
    chi_scf_err,
    chi_matrix_nscf,
    chi_nscf_err,
    chi_scf_inv_jacobs,
    chi_nscf_inv_jacobs,
):
    """
    Function to compute Hubbard U and Hund J values using scaled 2x2 formula,
    in addition to the associated uncertainty values
    - based on the study by Linscott et. al.
    """

    nx2 = 2 * site_index
    np1x2 = nx2 + 2

    # helpful functions for computing derivatives of "f" (Dyson) matrix
    def fmat_deriv_scf(kk, ll, ik, il):
        fd = chi_scf_inv_jacobs[nx2 + kk][nx2 + ll][nx2 + ik, nx2 + il]
        return fd

    def fmat_deriv_nscf(kk, ll, ik, il):
        fd = -chi_nscf_inv_jacobs[nx2 + kk][nx2 + ll][nx2 + ik, nx2 + il]
        return fd

    fmat = f_matrix[nx2:np1x2, nx2:np1x2]

    chi_sub_scf = chi_matrix_scf[nx2:np1x2, nx2:np1x2]
    chi_sub_scf_err = chi_scf_err[nx2:np1x2, nx2:np1x2]
    chi_sub_nscf_err = chi_nscf_err[nx2:np1x2, nx2:np1x2]

    # compute U value and error
    lam = (chi_sub_scf[0, 0] + chi_sub_scf[0, 1]) / (
        chi_sub_scf[1, 0] + chi_sub_scf[1, 1]
    )
    lam_deriv = np.zeros([2, 2])
    lam_deriv[0, 0] = 1.0 / (chi_sub_scf[1, 0] + chi_sub_scf[1, 1])
    lam_deriv[0, 1] = lam_deriv[0, 0]
    lam_deriv[1, 1] = -lam / (chi_sub_scf[1, 0] + chi_sub_scf[1, 1])
    lam_deriv[1, 0] = lam_deriv[1, 1]

    uval = (
        0.5 * (lam * (fmat[0, 0] + fmat[1, 0]) + fmat[0, 1] + fmat[1, 1]) / (lam + 1.0)
    )

    uval_err = 0.0
    # scf component
    u_deriv_scf = np.zeros([2, 2])
    for ik in [0, 1]:
        for il in [0, 1]:
            u_deriv_scf[ik, il] = (
                0.5
                / (lam + 1.0)
                * (
                    lam_deriv[ik, il] * (fmat[0, 0] + fmat[1, 0])
                    + lam
                    * (fmat_deriv_scf(0, 0, ik, il) + fmat_deriv_scf(1, 0, ik, il))
                    + fmat_deriv_scf(0, 1, ik, il)
                    + fmat_deriv_scf(1, 1, ik, il)
                )
                - uval / (lam + 1.0) * lam_deriv[ik, il]
            )
    jacob_vec = np.reshape(u_deriv_scf, [4, 1])
    uval_err = uval_err + np.sum(
        np.dot(
            np.transpose(jacob_vec),
            np.dot(np.diag(np.reshape(chi_sub_scf_err ** 2, [4])), jacob_vec),
        )
    )
    # nscf component
    u_deriv_nscf = np.zeros([2, 2])
    for ik in [0, 1]:
        for il in [0, 1]:
            u_deriv_nscf[ik, il] = (
                0.5
                / (lam + 1.0)
                * (
                    +lam
                    * (fmat_deriv_nscf(0, 0, ik, il) + fmat_deriv_nscf(1, 0, ik, il))
                    + fmat_deriv_nscf(0, 1, ik, il)
                    + fmat_deriv_nscf(1, 1, ik, il)
                )
            )
    jacob_vec = np.reshape(u_deriv_nscf, [4, 1])
    uval_err = uval_err + np.sum(
        np.dot(
            np.transpose(jacob_vec),
            np.dot(np.diag(np.reshape(chi_sub_nscf_err ** 2, [4])), jacob_vec),
        )
    )
    # compute std
    uval_err = np.sqrt(uval_err)

    # compute J value and error
    lam = (chi_sub_scf[0, 0] - chi_sub_scf[0, 1]) / (
        chi_sub_scf[1, 0] - chi_sub_scf[1, 1]
    )
    lam_deriv = np.zeros([2, 2])
    lam_deriv[0, 0] = 1.0 / (chi_sub_scf[1, 0] - chi_sub_scf[1, 1])
    lam_deriv[0, 1] = -lam_deriv[0, 0]
    lam_deriv[1, 1] = lam / (chi_sub_scf[1, 0] - chi_sub_scf[1, 1])
    lam_deriv[1, 0] = -lam_deriv[1, 1]

    jval = (
        -0.5 * (lam * (fmat[0, 0] - fmat[1, 0]) + fmat[0, 1] - fmat[1, 1]) / (lam - 1.0)
    )

    jval_err = 0.0
    # scf component
    j_deriv_scf = np.zeros([2, 2])
    for ik in [0, 1]:
        for il in [0, 1]:
            u_deriv_scf[ik, il] = (
                -0.5
                / (lam - 1.0)
                * (
                    lam_deriv[ik, il] * (fmat[0, 0] - fmat[1, 0])
                    + lam
                    * (fmat_deriv_scf(0, 0, ik, il) - fmat_deriv_scf(1, 0, ik, il))
                    + fmat_deriv_scf(0, 1, ik, il)
                    - fmat_deriv_scf(1, 1, ik, il)
                )
                + uval / (lam - 1.0) * lam_deriv[ik, il]
            )
    jacob_vec = np.reshape(j_deriv_scf, [4, 1])
    jval_err = jval_err + np.sum(
        np.dot(
            np.transpose(jacob_vec),
            np.dot(np.diag(np.reshape(chi_sub_scf_err ** 2, [4])), jacob_vec),
        )
    )
    # nscf component
    j_deriv_nscf = np.zeros([2, 2])
    for ik in [0, 1]:
        for il in [0, 1]:
            j_deriv_nscf[ik, il] = (
                -0.5
                / (lam - 1.0)
                * (
                    +lam
                    * (fmat_deriv_nscf(0, 0, ik, il) - fmat_deriv_nscf(1, 0, ik, il))
                    + fmat_deriv_nscf(0, 1, ik, il)
                    - fmat_deriv_nscf(1, 1, ik, il)
                )
            )
    jacob_vec = np.reshape(j_deriv_nscf, [4, 1])
    jval_err = jval_err + np.sum(
        np.dot(
            np.transpose(jacob_vec),
            np.dot(np.diag(np.reshape(chi_sub_nscf_err ** 2, [4])), jacob_vec),
        )
    )
    # compute std
    jval_err = np.sqrt(jval_err)

    return uval, uval_err, jval, jval_err
