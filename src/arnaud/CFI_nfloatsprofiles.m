

disp(['2022: ',num2str((4*diff(datenum(2022,[4,13],1))+3*diff(datenum(2022,[5,13],1))+3*diff(datenum(2022,[9,13],1)))/5,'%1.0f')])

disp(['2023: ',num2str((10*diff(datenum(2023,[1,13],1))+3*diff(datenum(2023,[4,13],1))+4*diff(datenum(2023,[5,13],1))+3*diff(datenum(2023,[9,13],1)))/5,'%1.0f')])

disp(['2024: ',num2str((20*diff(datenum(2024,[1,13],1))+3*diff(datenum(2024,[4,13],1))+4*diff(datenum(2024,[5,13],1))+3*diff(datenum(2024,[9,13],1)))/5,'%1.0f')])

disp(['2025: ',num2str(30*diff(datenum(2025,[1,13],1))/5,'%1.0f')])
