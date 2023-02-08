'''Sub-module to perform HDR tonemapping according to Fattal et al. 2002

https://www.cs.huji.ac.il/~danix/hdr/hdrc.pdf
'''
import cv2
import numpy as np
import scipy


def create_pyramids(imgdata, levels):
    pyrs = [imgdata]
    pyrdata = imgdata.copy()
    for i in range(levels):
        pyrdata = cv2.pyrDown(pyrdata)
        if (pyrdata.shape[0] >= 32) and (pyrdata.shape[1] >= 32):
            pyrs.append(pyrdata)
        else:
            print('Coarsest level of 32 reached. Stopping.')
            break
    return pyrs


def gradient(data, k):
    ''' Calculate the gradient using the central difference method. This is numpy.gradient.'''
    gradx = np.ones(data.shape)
    grady = np.ones(data.shape)
    for i in range(data.shape[1]):
        for j in range(data.shape[0]):
            if j-1 < 0:
                grady[i, j] = (data[i, j+1] - data[i, j]) / 2**k
            elif j+1 < data.shape[1]:
                grady[i, j] = (data[i, j+1] - data[i, j-1]) / 2**(k+1)
            else:
                grady[i, j] = (data[i, j] - data[i, j-1]) / 2**k

            if i-1 < 0:
                gradx[i, j] = (data[i+1, j] - data[i, j]) / 2**k
            elif i+1 < data.shape[0]:
                gradx[i, j] = (data[i+1, j] - data[i-1, j]) / 2**(k+1)
            else:
                gradx[i, j] = (data[i, j] - data[i-1, j]) / 2**k
    return [gradx, grady]


def phi(alpha, beta, gradnorm):
    ''' Scaling factor based on the magnitude of the gradient.'''
    return (alpha / gradnorm) * (gradnorm / alpha) ** beta


def gradient_attenuation(data, alpha, beta, levels=6):
    pyramids = create_pyramids(data, levels)
    # print(len(pyramids))
    grads = []
    phis = []
    # Initialise the cumulative attenuation factor to the same shape as the original image.
    atfac_total = np.ones_like(pyramids[0])
    for k, p in enumerate(pyramids):
        # print('Data shape:')
        # print(data.shape)
        # print('Pyramid shape:')
        # print(p.shape)
        grad = np.gradient(p, 2**k)
        grads.append(grad)
        print('Gradient is')
        print(grad)
        print('Gradient average is:')
        print(np.average(grad))
        print('Gradient norm is:')
        print(np.linalg.norm(grad, axis=0, ord=2))
        atfac = phi(0.1 * np.abs(np.average(grad)), beta, np.linalg.norm(grad, axis=0, ord=2))
        phis.append(atfac)
        print('Attenuation factor is')
        print(atfac)
        # print('Repeats is', data.shape[0] / p.shape[0])
        # print(np.repeat(atfac, repeats=data.shape[0] / p.shape[0], axis=0).repeat(data.shape[1] / p.shape[1], axis=1).shape)
        #atfac_total *= np.repeat(atfac, repeats=data.shape[0] / p.shape[0], axis=0).repeat(data.shape[1] / p.shape[1], axis=1)
        print('Interpolated phi:')
        print(atfac.shape)
        print(atfac_total.shape)
        print(cv2.resize(atfac, dsize=atfac_total.shape[::-1], interpolation=cv2.INTER_LINEAR).shape)
        atfac_total *= cv2.resize(atfac, dsize=atfac_total.shape[::-1], interpolation=cv2.INTER_LINEAR)
        return atfac_total
    return atfac_total


def gradH(data):
    ''' Calculate the gradient using the central difference method. This is numpy.gradient.'''
    gradx = np.ones(data.shape)
    grady = np.ones(data.shape)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if j+1 < data.shape[1]:
                grady[i, j] = (data[i, j+1] - data[i, j])
            else:
                grady[i, j] = (data[i, j] - data[i, j])

            if i+1 < data.shape[0]:
                gradx[i, j] = (data[i+1, j] - data[i, j])
            else:
                gradx[i, j] = (data[i, j] - data[i, j])
    return [gradx, grady]


def divG(data):
    ''' Calculate the divergence of G using backward difference.'''
    div = np.ones_like(data[0])
    for i in range(data[0].shape[0]):
        for j in range(data[0].shape[1]):
            div[i, j] = data[0][i, j] - data[0][i-1, j] + data[1][i, j] - data[1][i, j-1]
    return div


def map_fattal(data, alpha, beta, levels=6):
    plot_matrix(data, 'fattal_step1_data.png')
    # Calculate the gradient attenuation.
    attenuation = gradient_attenuation(data, alpha, beta, levels=levels)
    plot_matrix(attenuation, 'fattal_step2_attenuation.png')
    plot_matrix(attenuation*data, 'fattal_step3_attenuation_x_data.png')
    # print(attenuation)
    # Compress the gradient to reduce the dynamic range.
    data_grad = np.gradient(data)
    # print(data_grad)n
    # print(compgrad)
    data_grad = gradH(data)
    data_grad[0][0, :] = 0
    # data_grad[0][:, 0] = 0
    # data_grad[0][:, -1] = 0
    data_grad[0][-1, :] = 0
    
    # data_grad[1][0, :] = 0
    data_grad[1][:, 0] = 0
    data_grad[1][:, -1] = 0
    # data_grad[1][-1, :] = 0
    plot_matrix(data_grad[0], 'fattal_step4_gradX.png')
    plot_matrix(data_grad[1], 'fattal_step4_gradY.png')
    G = np.asarray(data_grad) * 0.1#attenuation
    G_div = divG(G)
    plot_matrix(G_div, 'fattal_step5_divG.png')
    # Finally recover the reduced dynamic range image.
    # Solve the Poisson equation
    # https://folk.ntnu.no/leifh/teaching/tkt4140/._main055.html
    # we'll solve Ax=b where x is the vector of I[i,j] and b is the right hand side
    N = data.shape[0] * data.shape[1]
    d = np.ones(N)
    diag = -4 * d
    diag1 = d[:-1]
    diag5 = d[:N-5]
    # The matrix is mostly empty, use a sparse matrix.
    A = scipy.sparse.diags([diag, diag1, diag1, diag5, diag5], offsets=[0, 1, -1, 5, -5], format='csc')

    # The vector "b".
    print('First row of div G:')
    print(G_div[0, :])
    b_temp = G_div.copy()
    # Set the boundary conditions
    b_temp[0, :] = 0
    b_temp[:, 0] = 0
    b_temp[:, -1] = 0
    b_temp[-1, :] = 0
    plot_matrix(b_temp, 'fattal_step5_divG_b.png')
    b = np.matrix(b_temp).flatten().T
    print('Matrix:')
    # # print(A.toarray())
    # plot_matrix(A.toarray(), 'fattal_step6_matrix.png')
    I = scipy.sparse.linalg.spsolve(A, b).reshape(data.shape)
    plot_matrix(I, 'fattal_step6_I.png')
    return np.exp(I)

from matplotlib.pyplot import figure
def plot_matrix(data, name):
    fig = figure()
    ax = fig.add_subplot(111)
    ax.imshow(data)
    fig.savefig(name, dpi=300)

np.random.seed(1231292921)
N = 1024
x = np.random.rand(N, N)
# from astropy.io import fits
# x = fits.getdata('tests_output/LoTSS-DR2_202.4842_47.2306_0.300.fits').squeeze()
# y = cv2.imread('/home/frits/EveryStamp/tests_output/LoTSS-DR2_202.4842_47.2306_0.300_linear.png')
# x = cv2.cvtColor(y, cv2.COLOR_BGR2GRAY).astype(float)
print(x.shape)
map_fattal(x, 0.1, 0.85)
# x2 = cv2.pyrDown(x)
# print(x2.shape)
# print(cv2.pyrUp(x2).shape)
# y = np.linalg.norm(np.gradient(x), axis=0)
# print(y)
# print(phi(0.1 * np.average(y), 0.85, y))

# print(x.shape)
# print(x2.shape)
# x3 = cv2.pyrDown(x2)
# print(gradient(x, 0)[0] - np.gradient(x)[0], '\n', gradient(x, 0)[1] - np.gradient(x)[1])
# print(gradient(x2, 1)[0] - np.gradient(x2, 2)[0], '\n', gradient(x2, 1)[1] - np.gradient(x2, 2)[1])
# print(gradient(x3, 2)[0] - np.gradient(x3, 4)[0], '\n', gradient(x3, 2)[1] - np.gradient(x3, 4)[1])