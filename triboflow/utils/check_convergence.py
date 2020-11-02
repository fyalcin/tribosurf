
class CheckConvergence():
    
    """
    This class contains the methods for checking the convergence.
    """

    def IsListConverged(input_list, tol, n=3):
        """Check if the last n values of an array are within tol of each other.
        

        Parameters
        ----------
        input_list : list of float
            Total energies to be checked for convergence
        tol : float
            Tolerance for the convergence.
        n : int, optional
            Number of entries at the end of energy_list that have to be within
            etol for the list to be considered converged. The default is 3.

        Returns
        -------
        Bool
            True if input_list is converged, False otherwise.

        """
        if len(input_list) <= n:
            return False
        else:
            check_list = [False]*n
            l = input_list.copy()
            l.reverse()
            for i, b in enumerate(check_list):
                if abs(l[0]-l[i+1]) < tol:
                    check_list[i] = True
            return all(check_list)