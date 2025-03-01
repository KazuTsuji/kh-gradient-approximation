function convertFigToPDF(figName, outPutPDFName)
    % .figファイルを開く
    f = openfig(figName, 'invisible'); % 'invisible'で非表示のまま開く
    
    % 出力サイズをfigureサイズに合わせる
    f.Units = 'centimeters';
    f.PaperUnits = f.Units;
    f.PaperPosition = [0, 0, f.Position(3:4)];
    f.PaperSize = f.Position(3:4);
    
    % figureをPDFで出力
    print(f, outPutPDFName, '-dpdf');
    
    % figureを閉じる
    close(f);
