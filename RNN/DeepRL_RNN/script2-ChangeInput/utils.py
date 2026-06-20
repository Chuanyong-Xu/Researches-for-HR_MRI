import random

import numpy as np
import torch


def set_seed(seed=42):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


#def reshape_and_mask_tdev(diff, switch):
#    """
#    先调用 reshape_own，对 tf_all 进行判断：
#    如果某个 trial 的 tf_all == [0, 1]，则将对应的 tdev_all 改为 [0, 0, 0]
#
#    :param switch:
#    :param diff:
#    :return: 处理后的 tdev_all, tf_all
#    """
#
#    # 创建 mask，判断 tf_all 是否为 [0, 1]
#    mask = torch.all(switch == torch.tensor([0.0, 1.0], device=switch.device), dim=1)
#
#    # 将对应的 tdev_all 设置为 [0, 0, 0]
#    diff[mask] = torch.tensor([0.0, 0.0, 0.0], device=diff.device)
#
#    return diff, switch


def reshape_and_mask_tdev(outcome, switch):
    """
    outcome: (B, T, 6)
    switch : (B, T, 2)  one-hot, [0,1] or [1,0]
    当 switch == [0,1] 时，将 outcome 对应位置置零
    """
    # mask: (B, T)  True 表示 tf == [0,1]
    mask = (switch[..., 0] == 0) & (switch[..., 1] == 1)  # (B,T)

    # outcome[mask] 会变成 (N_true, 6)，可直接置 0
    outcome = outcome.clone()  # 更安全，避免原地改导致 autograd 问题
    outcome[mask] = 0.0
    return outcome