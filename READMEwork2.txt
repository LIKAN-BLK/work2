��� ����� ��� �����: https://mega.co.nz/#!uY8mATYI!VkwHVxaEQoIbh2P-b2EQ7TcMGUVh3-hjmNgRn06qzhE
���� ����� ������ ��� ��������� �����  (06 � 07)

�������� ������� http://bcilab.atedgeof.net/wiki/%D0%9C%D0%B5%D1%82%D0%BA%D0%B8_%D0%B0%D0%B9%D0%BB%D0%B0%D0%B9%D0%BD%D0%B7


>> addpath('C:\Users\LIKAN\BCI\work2')
>> vars_e206

>> [eegT, eegNT, parameters, fixationDuration, sRate, path, epoch_size] = eye_loaddata('C:\Users\LIKAN\BCI\e2_all_files\06', 400, 400, 0);

path - ����, ��� ����� �����
fixation_threshold - ���������� ����� �� ���� � ���� ������� ��������
epoch_size - ������ �����
left_border - ����� ������� ����� ������ ������� �������� ����� ���


>> params = eye_train(eegT, eegNT, chans_e206, parameters, fixationDuration, sRate, path, epoch_size);


������ eeg (����� ������� �� �������, �����, ����� ������������ (?))